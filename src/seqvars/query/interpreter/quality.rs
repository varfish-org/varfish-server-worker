use crate::{
    common::strip_gt_leading_slash,
    seqvars::query::schema::{CallInfo, CaseQuery, FailChoice, QualitySettings, SequenceVariant},
};

/// Return type for the `passes` function.
#[derive(Debug)]
pub struct PassOrNoCall {
    /// Whether the variant should be kept.
    pub pass: bool,
    /// For which samples should the genotype be interpreted as no-call.
    pub no_call_samples: Vec<String>,
}

/// Determine whether the `SequenceVariant` passes the quality filter.
/// Will return `FailChoice::Ignore` if the variant passes.
pub fn passes(query: &CaseQuery, seqvar: &SequenceVariant) -> Result<PassOrNoCall, anyhow::Error> {
    let mut result = PassOrNoCall {
        pass: true,
        no_call_samples: Vec::new(),
    };
    for (sample_name, quality_settings) in &query.quality {
        if let Some(call_info) = seqvar.call_info.get(sample_name) {
            if let Some(fail) = passes_for_sample(quality_settings, call_info) {
                match fail {
                    FailChoice::Ignore => {
                        // ignore quality failure for sample
                    }
                    FailChoice::Drop => {
                        tracing::trace!(
                            "sample {} (call_info={:?}) in variant {:?} fails quality filter {:?}",
                            &sample_name,
                            &call_info,
                            &seqvar,
                            &quality_settings
                        );
                        result.pass = false;
                        break;
                    }
                    FailChoice::NoCall => {
                        result.no_call_samples.push(sample_name.clone());
                    }
                }
            } else {
                // no failure, all good
            }
        } else {
            anyhow::bail!("sample {} not found in call info", sample_name);
        }
    }

    tracing::trace!(
        "variant {:?} has result {:?} for quality filter {:?}",
        seqvar,
        result,
        &query.genotype
    );
    Ok(result)
}

/// Return failure code (or None for all-pass) for one sample's call info.
fn passes_for_sample(
    quality_settings: &QualitySettings,
    call_info: &CallInfo,
) -> Option<FailChoice> {
    #[derive(PartialEq, Eq)]
    enum Genotype {
        Het,
        Hom,
        Ref,
        NoCall,
    }

    let genotype = if let Some(genotype) = call_info.genotype.as_ref() {
        let genotype = strip_gt_leading_slash(genotype);
        match genotype {
            "0/1" | "1/0" | "0|1" | "1|0" => Genotype::Het,
            "1/1" | "1|1" | "1" => Genotype::Hom,
            "0/0" | "0|0" | "0" => Genotype::Ref,
            _ => Genotype::NoCall,
        }
    } else {
        Genotype::NoCall
    };

    // dp_het/dp_hom and ab
    match genotype {
        Genotype::Het => {
            // dp_het
            if let Some(dp_het) = quality_settings.dp_het {
                if let Some(dp) = call_info.dp {
                    if dp < dp_het {
                        return Some(quality_settings.fail);
                    }
                }
            }

            // ab
            if let (Some(settings_ab), Some(call_dp), Some(call_ad)) =
                (quality_settings.ab, call_info.dp, call_info.ad)
            {
                let ab_raw = call_ad as f64 / call_dp as f64;
                let ab = if ab_raw > 0.5 { 1.0 - ab_raw } else { ab_raw };
                let eps = 1e-6f64;
                if ab + eps < settings_ab as f64 {
                    return Some(quality_settings.fail);
                }
            }
        }
        Genotype::Hom => {
            if let Some(dp_hom) = quality_settings.dp_hom {
                if let Some(dp) = call_info.dp {
                    if dp < dp_hom {
                        return Some(quality_settings.fail);
                    }
                }
            }
        }
        Genotype::Ref | Genotype::NoCall => (),
    }

    // gq
    if let (Some(settings_gq), Some(call_gq)) = (quality_settings.gq, call_info.quality) {
        if call_gq < settings_gq as f32 {
            return Some(quality_settings.fail);
        }
    }

    if genotype != Genotype::Ref {
        // ad
        if let (Some(settings_ad), Some(call_ad)) = (quality_settings.ad, call_info.ad) {
            if call_ad < settings_ad {
                return Some(quality_settings.fail);
            }
        }

        // ad_max
        if let (Some(settings_ad_max), Some(call_ad)) = (quality_settings.ad_max, call_info.ad) {
            if call_ad > settings_ad_max {
                return Some(quality_settings.fail);
            }
        }
    }

    None
}

#[cfg(test)]
mod test {
    use rstest::rstest;

    use crate::seqvars::query::schema::{
        CallInfo, CaseQuery,
        FailChoice::{self, *},
        QualitySettings, SequenceVariant,
    };

    #[rstest]
    #[case(Ignore, true, true, false)]
    #[case(Ignore, false, true, false)]
    #[case(Drop, true, true, false)]
    #[case(Drop, false, false, false)]
    #[case(NoCall, true, true, false)]
    #[case(NoCall, false, true, true)]
    fn passes(
        #[case] q_fail: FailChoice,
        #[case] should_pass: bool,
        #[case] expected_pass: bool,
        #[case] any_no_call_sample: bool,
    ) -> Result<(), anyhow::Error> {
        let query = CaseQuery {
            quality: vec![(
                String::from("sample"),
                QualitySettings {
                    dp_het: None,
                    dp_hom: None,
                    gq: if should_pass { None } else { Some(40) },
                    ab: None,
                    ad: None,
                    ad_max: None,
                    fail: q_fail,
                },
            )]
            .into_iter()
            .collect(),
            ..Default::default()
        };
        let seqvar = SequenceVariant {
            call_info: vec![(
                String::from("sample"),
                CallInfo {
                    quality: if should_pass { None } else { Some(30f32) },
                    ..Default::default()
                },
            )]
            .into_iter()
            .collect(),
            ..Default::default()
        };

        let res = super::passes(&query, &seqvar)?;

        assert_eq!(
            res.pass, expected_pass,
            "query = {:#?}, seqvar = {:#?}",
            &query, &seqvar
        );
        let expected = if any_no_call_sample {
            vec![String::from("sample")]
        } else {
            vec![]
        };
        assert_eq!(
            &res.no_call_samples, &expected,
            "query = {:#?}, seqvar = {:#?}",
            &query, &seqvar
        );

        Ok(())
    }

    #[rstest]
    #[allow(clippy::too_many_arguments)]
    // het, pass dp
    #[case(
        Some(10), // q_dp_het
        None, // q_dp_hom
        None, // q_gq
        None, // q_ab
        None, // q_ad
        None, // q_ad_max
        Ignore, // q_fail
        Some("0/1"), // c_genotype
        None, // c_quality
        Some(10), // c_dp
        None, // c_ad
        None,  // expected, None = pass
    )]
    // het, fail dp
    #[case(
        Some(10), // q_dp_het
        None, // q_dp_hom
        None, // q_gq
        None, // q_ab
        None, // q_ad
        None, // q_ad_max
        Ignore, // q_fail
        Some("0/1"), // c_genotype
        None, // c_quality
        Some(9), // c_dp
        None, // c_ad
        Some(Ignore),  // expected, None = pass
    )]
    // hom, pass dp
    #[case(
        None, // q_dp_het
        Some(10), // q_dp_hom
        None, // q_gq
        None, // q_ab
        None, // q_ad
        None, // q_ad_max
        Ignore, // q_fail
        Some("1/1"), // c_genotype
        None, // c_quality
        Some(10), // c_dp
        None, // c_ad
        None,  // expected, None = pass
    )]
    // hom, fail dp
    #[case(
        None, // q_dp_het
        Some(10), // q_dp_hom
        None, // q_gq
        None, // q_ab
        None, // q_ad
        None, // q_ad_max
        Ignore, // q_fail
        Some("1/1"), // c_genotype
        None, // c_quality
        Some(9), // c_dp
        None, // c_ad
        Some(Ignore),  // expected, None = pass
    )]
    // pass gq
    #[case(
        None, // q_dp_het
        None, // q_dp_hom
        Some(10), // q_gq
        None, // q_ab
        None, // q_ad
        None, // q_ad_max
        Ignore, // q_fail
        None, // c_genotype
        Some(10f32), // c_quality
        None, // c_dp
        None, // c_ad
        None,  // expected, None = pass
    )]
    // fail gq
    #[case(
        None, // q_dp_het
        None, // q_dp_hom
        Some(10), // q_gq
        None, // q_ab
        None, // q_ad
        None, // q_ad_max
        Ignore, // q_fail
        Some("0/1"), // c_genotype
        Some(9f32), // c_quality
        None, // c_dp
        None, // c_ad
        Some(Ignore),  // expected, None = pass
    )]
    // het, pass ab lower
    #[case(
        None, // q_dp_het
        None, // q_dp_hom
        None, // q_gq
        Some(0.2), // q_ab
        None, // q_ad
        None, // q_ad_max
        Ignore, // q_fail
        Some("0/1"), // c_genotype
        None, // c_quality
        Some(100), // c_dp
        Some(20), // c_ad
        None,  // expected, None = pass
    )]
    // het, pass ab upper
    #[case(
        None, // q_dp_het
        None, // q_dp_hom
        None, // q_gq
        Some(0.2), // q_ab
        None, // q_ad
        None, // q_ad_max
        Ignore, // q_fail
        Some("0/1"), // c_genotype
        None, // c_quality
        Some(100), // c_dp
        Some(80), // c_ad
        None,  // expected, None = pass
    )]
    // het, fail ab lower
    #[case(
        None, // q_dp_het
        None, // q_dp_hom
        None, // q_gq
        Some(0.2), // q_ab
        None, // q_ad
        None, // q_ad_max
        Ignore, // q_fail
        Some("0/1"), // c_genotype
        None, // c_quality
        Some(100), // c_dp
        Some(19), // c_ad
        Some(Ignore),  // expected, None = pass
    )]
    // het, fail ab upper
    #[case(
        None, // q_dp_het
        None, // q_dp_hom
        None, // q_gq
        Some(0.2), // q_ab
        None, // q_ad
        None, // q_ad_max
        Ignore, // q_fail
        Some("0/1"), // c_genotype
        None, // c_quality
        Some(100), // c_dp
        Some(81), // c_ad
        Some(Ignore),  // expected, None = pass
    )]
    // hom, ab ignored
    #[case(
        None, // q_dp_het
        None, // q_dp_hom
        None, // q_gq
        Some(0.2), // q_ab
        None, // q_ad
        None, // q_ad_max
        Ignore, // q_fail
        Some("1/1"), // c_genotype
        None, // c_quality
        None, // c_dp
        None, // c_ad
        None,  // expected, None = pass
    )]
    // pass ad
    #[case(
        None, // q_dp_het
        None, // q_dp_hom
        None, // q_gq
        None, // q_ab
        Some(10), // q_ad
        None, // q_ad_max
        Ignore, // q_fail
        None, // c_genotype
        None, // c_quality
        None, // c_dp
        Some(10), // c_ad
        None,  // expected, None = pass
    )]
    // fail ad
    #[case(
        None, // q_dp_het
        None, // q_dp_hom
        None, // q_gq
        None, // q_ab
        Some(10), // q_ad
        None, // q_ad_max
        Ignore, // q_fail
        None, // c_genotype
        None, // c_quality
        None, // c_dp
        Some(9), // c_ad
        Some(Ignore),  // expected, None = pass
    )]
    // pass ad_max
    #[case(
        None, // q_dp_het
        None, // q_dp_hom
        None, // q_gq
        None, // q_ab
        None, // q_ad
        Some(10), // q_ad_max
        Ignore, // q_fail
        None, // c_genotype
        None, // c_quality
        None, // c_dp
        Some(10), // c_ad
        None,  // expected, None = pass
    )]
    // fail ad_max
    #[case(
        None, // q_dp_het
        None, // q_dp_hom
        None, // q_gq
        None, // q_ab
        None, // q_ad
        Some(10), // q_ad_max
        Ignore, // q_fail
        None, // c_genotype
        None, // c_quality
        None, // c_dp
        Some(11), // c_ad
        Some(Ignore),  // expected, None = pass
    )]
    // all none
    #[case(
        None, // q_dp_het
        None, // q_dp_hom
        None, // q_gq
        None, // q_ab
        None, // q_ad
        None, // q_ad_max
        Ignore, // q_fail
        None, // c_genotype
        None, // c_quality
        None, // c_dp
        None, // c_ad
        None,  // expected, None = pass
    )]
    fn passes_for_sample(
        #[case] q_dp_het: Option<i32>,
        #[case] q_dp_hom: Option<i32>,
        #[case] q_gq: Option<i32>,
        #[case] q_ab: Option<f32>,
        #[case] q_ad: Option<i32>,
        #[case] q_ad_max: Option<i32>,
        #[case] q_fail: FailChoice,
        #[case] c_genotype: Option<&'static str>,
        #[case] c_quality: Option<f32>,
        #[case] c_dp: Option<i32>,
        #[case] c_ad: Option<i32>,
        #[case] expected: Option<FailChoice>,
    ) -> Result<(), anyhow::Error> {
        let settings = QualitySettings {
            dp_het: q_dp_het,
            dp_hom: q_dp_hom,
            gq: q_gq,
            ab: q_ab,
            ad: q_ad,
            ad_max: q_ad_max,
            fail: q_fail,
        };
        let call_info = CallInfo {
            genotype: c_genotype.map(|s| s.to_string()),
            quality: c_quality,
            dp: c_dp,
            ad: c_ad,
            ..Default::default()
        };

        assert_eq!(
            super::passes_for_sample(&settings, &call_info),
            expected,
            "settings: {:?}, call info: {:?}",
            settings,
            call_info
        );

        Ok(())
    }
}
