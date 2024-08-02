use crate::{
    common::strip_gt_leading_slash,
    seqvars::query::schema::{
        data::{CallInfo, VariantRecord},
        query::{CaseQuery, SampleQualitySettings},
    },
};

/// Return type for the `passes` function.
#[derive(Debug)]
pub struct PassOrNoCall {
    /// Whether the variant should be kept.
    pub pass: bool,
    /// For which samples should the genotype be interpreted as no-call.
    pub no_call_samples: Vec<String>,
}

/// Determine whether the `VariantRecord` passes the quality filter.
/// Will return `FailFilterChoice::Ignore` if the variant passes.
pub fn passes(query: &CaseQuery, seqvar: &VariantRecord) -> Result<PassOrNoCall, anyhow::Error> {
    let mut result = PassOrNoCall {
        pass: true,
        no_call_samples: Vec::new(),
    };
    for (sample_name, quality_settings) in &query.quality.sample_qualities {
        if let Some(call_info) = seqvar.call_infos.get(sample_name) {
            if !passes_for_sample(quality_settings, call_info) {
                tracing::trace!(
                    "sample {} (call_info={:?}) in variant {:?} fails quality filter {:?}",
                    &sample_name,
                    &call_info,
                    &seqvar,
                    &quality_settings
                );
                if quality_settings.filter_active {
                    result.pass = false;
                    break;
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

/// Return whether the sample passes the quality filter.
fn passes_for_sample(quality_settings: &SampleQualitySettings, call_info: &CallInfo) -> bool {
    // Ad-hoc enum for genotype.
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

    // min_dp_het, min_dp_hom, and min_ab
    match genotype {
        Genotype::Het => {
            // min_dp_het
            if let (Some(min_dp_het), Some(dp)) = (quality_settings.min_dp_het, call_info.dp) {
                if dp < min_dp_het {
                    return false;
                }
            }

            // min_ab
            if let (Some(min_ab), Some(call_dp), Some(call_ad)) =
                (quality_settings.min_ab, call_info.dp, call_info.ad)
            {
                let ab_raw = call_ad as f64 / call_dp as f64;
                let ab = if ab_raw > 0.5 { 1.0 - ab_raw } else { ab_raw };
                let eps = 1e-6f64;
                if ab + eps < min_ab as f64 {
                    return false;
                }
            }
        }
        Genotype::Hom => {
            // min_dp_hom
            if let (Some(dp_hom), Some(dp)) = (quality_settings.min_dp_hom, call_info.dp) {
                if dp < dp_hom {
                    return false;
                }
            }
        }
        Genotype::Ref | Genotype::NoCall => (),
    }

    // min_gq
    if let (Some(settings_gq), Some(call_gq)) = (quality_settings.min_gq, call_info.quality) {
        if call_gq < settings_gq as f32 {
            return false;
        }
    }

    if genotype != Genotype::Ref {
        // min_ad
        if let (Some(settings_ad), Some(call_ad)) = (quality_settings.min_ad, call_info.ad) {
            if call_ad < settings_ad {
                return false;
            }
        }

        // max_ad
        if let (Some(settings_ad_max), Some(call_ad)) = (quality_settings.max_ad, call_info.ad) {
            if call_ad > settings_ad_max {
                return false;
            }
        }
    }

    true
}

#[cfg(test)]
mod test {
    use crate::seqvars::query::schema::data::{CallInfo, VariantRecord};
    use crate::seqvars::query::schema::query::{
        CaseQuery, QuerySettingsQuality, SampleQualitySettings,
    };

    #[rstest::rstest]
    #[case(false, true, true, false)]
    #[case(false, false, true, false)]
    #[case(true, true, true, false)]
    #[case(true, false, false, false)]
    fn passes(
        #[case] filter_active: bool,
        #[case] should_pass: bool,
        #[case] expected_pass: bool,
        #[case] any_no_call_sample: bool,
    ) -> Result<(), anyhow::Error> {
        let query = CaseQuery {
            quality: QuerySettingsQuality {
                sample_qualities: indexmap::indexmap! {
                    String::from("sample") =>
                    SampleQualitySettings {
                        sample: String::from("sample"),
                        filter_active,
                        min_gq: if should_pass { None } else { Some(40) },
                        ..Default::default()
                    },
                },
            },
            ..Default::default()
        };
        let seqvar = VariantRecord {
            call_infos: indexmap::indexmap! {
                String::from("sample") =>
                CallInfo {
                    quality: if should_pass { None } else { Some(30f32) },
                    ..Default::default()
                },
            },
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

    #[rstest::rstest]
    #[allow(clippy::too_many_arguments)]
    // het, pass dp
    #[case(
        Some(10), // q_min_dp_het
        None, // q_dpmin__hom
        None, // q_min_gq
        None, // q_min_ab
        None, // q_min_ad
        None, // q_max_ad
        false, // filter_active
        Some("0/1"), // c_genotype
        None, // c_quality
        Some(10), // c_dp
        None, // c_ad
        true, // expected
    )]
    // het, fail dp
    #[case(
        Some(10), // q_min_dp_het
        None, // q_dpmin__hom
        None, // q_min_gq
        None, // q_min_ab
        None, // q_min_ad
        None, // q_max_ad
        false, // filter_active
        Some("0/1"), // c_genotype
        None, // c_quality
        Some(9), // c_dp
        None, // c_ad
        true,  // expected
    )]
    // hom, pass dp
    #[case(
        None, // q_min_dp_het
        Some(10), // q_min_dp_hom
        None, // q_min_gq
        None, // q_min_ab
        None, // q_min_ad
        None, // q_max_ad
        false, // filter_active
        Some("1/1"), // c_genotype
        None, // c_quality
        Some(10), // c_dp
        None, // c_ad
        true, // expected
    )]
    // hom, fail dp
    #[case(
        None, // q_min_dp_het
        Some(10), // q_min_dp_hom
        None, // q_min_gq
        None, // q_min_ab
        None, // q_min_ad
        None, // q_max_ad
        false, // filter_active
        Some("1/1"), // c_genotype
        None, // c_quality
        Some(9), // c_dp
        None, // c_ad
        true,  // expected
    )]
    // pass gq
    #[case(
        None, // q_min_dp_het
        None, // q_min_dp_hom
        Some(10), // min_q_gq
        None, // q_min_ab
        None, // q_min_ad
        None, // q_max_ad
        false, // filter_active
        None, // c_genotype
        Some(10f32), // c_quality
        None, // c_dp
        None, // c_ad
        true, // expected
    )]
    // fail gq
    #[case(
        None, // q_min_dp_het
        None, // q_min_dp_hom
        Some(10), // min_q_gq
        None, // q_min_ab
        None, // q_min_ad
        None, // q_max_ad
        false, // filter_active
        Some("0/1"), // c_genotype
        Some(9f32), // c_quality
        None, // c_dp
        None, // c_ad
        true,  // expected
    )]
    // het, pass ab lower
    #[case(
        None, // q_min_dp_het
        None, // q_min_dp_hom
        None, // q_min_gq
        Some(0.2), //min_ q_ab
        None, // q_min_ad
        None, // q_max_ad
        false, // filter_active
        Some("0/1"), // c_genotype
        None, // c_quality
        Some(100), // c_dp
        Some(20), // c_ad
        true, // expected
    )]
    // het, pass ab upper
    #[case(
        None, // q_min_dp_het
        None, // q_min_dp_hom
        None, // q_min_gq
        Some(0.2), //min_ q_ab
        None, // q_min_ad
        None, // q_max_ad
        false, // filter_active
        Some("0/1"), // c_genotype
        None, // c_quality
        Some(100), // c_dp
        Some(80), // c_ad
        true, // expected
    )]
    // het, fail ab lower
    #[case(
        None, // q_min_dp_het
        None, // q_min_dp_hom
        None, // q_min_gq
        Some(0.2), //min_ q_ab
        None, // q_min_ad
        None, // q_max_ad
        false, // filter_active
        Some("0/1"), // c_genotype
        None, // c_quality
        Some(100), // c_dp
        Some(19), // c_ad
        true,  // expected
    )]
    // het, fail ab upper
    #[case(
        None, // q_min_dp_het
        None, // q_min_dp_hom
        None, // q_min_gq
        Some(0.2), //min_ q_ab
        None, // q_min_ad
        None, // q_max_ad
        false, // filter_active
        Some("0/1"), // c_genotype
        None, // c_quality
        Some(100), // c_dp
        Some(81), // c_ad
        true,  // expected
    )]
    // hom, ab ignored
    #[case(
        None, // q_min_dp_het
        None, // q_min_dp_hom
        None, // q_min_gq
        Some(0.2), //min_ q_ab
        None, // q_min_ad
        None, // q_max_ad
        false, // filter_active
        Some("1/1"), // c_genotype
        None, // c_quality
        None, // c_dp
        None, // c_ad
        true, // expected
    )]
    // pass ad
    #[case(
        None, // q_min_dp_het
        None, // q_min_dp_hom
        None, // q_min_gq
        None, // q_min_ab
        Some(10), // min_q_ad
        None, // q_max_ad
        false, // filter_active
        None, // c_genotype
        None, // c_quality
        None, // c_dp
        Some(10), // c_ad
        true, // expected
    )]
    // fail ad
    #[case(
        None, // q_min_dp_het
        None, // q_min_dp_hom
        None, // q_min_gq
        None, // q_min_ab
        Some(10), // min_q_ad
        None, // q_max_ad
        false, // filter_active
        None, // c_genotype
        None, // c_quality
        None, // c_dp
        Some(9), // c_ad
        true,  // expected
    )]
    // pass ad_max
    #[case(
        None, // q_min_dp_het
        None, // q_min_dp_hom
        None, // q_min_gq
        None, // q_min_ab
        None, // q_min_ad
        Some(10), // max_ad
        false, // filter_active
        None, // c_genotype
        None, // c_quality
        None, // c_dp
        Some(10), // c_ad
        true, // expected
    )]
    // fail ad_max
    #[case(
        None, // q_min_dp_het
        None, // q_min_dp_hom
        None, // q_min_gq
        None, // q_min_ab
        None, // q_min_ad
        Some(10), // max_ad
        false, // filter_active
        None, // c_genotype
        None, // c_quality
        None, // c_dp
        Some(11), // c_ad
        true,  // expected
    )]
    // all none
    #[case(
        None, // q_min_dp_het
        None, // q_min_dp_hom
        None, // q_min_gq
        None, // q_min_ab
        None, // q_min_ad
        None, // q_max_ad
        false, // filter_active
        None, // c_genotype
        None, // c_quality
        None, // c_dp
        None, // c_ad
        true, // expected
    )]
    fn passes_for_sample(
        #[case] q_min_dp_het: Option<i32>,
        #[case] q_min_dp_hom: Option<i32>,
        #[case] q_min_gq: Option<i32>,
        #[case] q_min_ab: Option<f32>,
        #[case] q_min_ad: Option<i32>,
        #[case] q_max_ad: Option<i32>,
        #[case] filter_active: bool,
        #[case] c_genotype: Option<&'static str>,
        #[case] c_quality: Option<f32>,
        #[case] c_dp: Option<i32>,
        #[case] c_ad: Option<i32>,
        #[case] expected: bool,
    ) -> Result<(), anyhow::Error> {
        let settings = SampleQualitySettings {
            sample: String::from("sample"),
            filter_active,
            min_dp_het: q_min_dp_het,
            min_dp_hom: q_min_dp_hom,
            min_gq: q_min_gq,
            min_ab: q_min_ab,
            min_ad: q_min_ad,
            max_ad: q_max_ad,
        };
        let call_info = CallInfo {
            sample: String::from("sample"),
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
