use crate::seqvars::query::schema::{
    CallInfo, CaseQuery, FailChoice, QualitySettings, SequenceVariant,
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
    enum Genotype {
        Het,
        Hom,
        Ref,
        NoCall,
    }

    let genotype = if let Some(genotype) = call_info.genotype.as_ref() {
        match genotype.as_str() {
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
                let ab_raw = call_ad as f32 / call_dp as f32;
                let ab = if ab_raw > 0.5 { 1.0 - ab_raw } else { ab_raw };
                if ab < settings_ab {
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

    None
}
