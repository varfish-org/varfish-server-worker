use crate::seqvars::query::schema::{CaseQuery, GenomicRegion, Range, SequenceVariant};

/// Determine whether the `SequenceVariant` passes the regions allowlist filter.
pub fn passes(query: &CaseQuery, seqvar: &SequenceVariant) -> bool {
    if let Some(region_allowlist) = &query.genomic_regions {
        if region_allowlist.is_empty() {
            true
        } else {
            let res = region_allowlist.iter().any(|region| {
                overlaps(
                    region,
                    &seqvar.chrom,
                    seqvar.pos,
                    seqvar.pos + seqvar.reference.len() as i32 - 1,
                )
            });
            if !res {
                tracing::trace!(
                    "variant {:?} fails region allowlist filter {:?}",
                    seqvar,
                    &region_allowlist
                );
            }
            res
        }
    } else {
        true
    }
}

fn overlaps(region: &GenomicRegion, seqvar_chrom: &str, seqvar_pos: i32, seqvar_end: i32) -> bool {
    let GenomicRegion {
        chrom: region_chrom,
        range: region_range,
    } = region;

    let region_chrom_c = annonars::common::cli::canonicalize(region_chrom);
    let seqvar_chrom_c = annonars::common::cli::canonicalize(seqvar_chrom);
    if region_chrom_c != seqvar_chrom_c {
        return false;
    }

    if let Some(Range {
        start: region_start,
        end: region_end,
    }) = region_range
    {
        *region_start <= seqvar_end && *region_end >= seqvar_pos
    } else {
        true
    }
}

#[cfg(test)]
mod test {
    use rstest::rstest;

    #[rstest]
    #[case("1", Some((100, 200)), "1", 100, 200, true)]
    #[case("chr1", Some((100, 200)), "1", 100, 200, true)]
    #[case("chr1", Some((100, 200)), "chr1", 100, 200, true)]
    #[case("chr1", Some((100, 200)), "chr1", 200, 300, true)]
    #[case("chr1", Some((100, 200)), "chr1", 201, 300, false)]
    #[case("1", Some((100, 200)), "2", 100, 200, false)]
    fn overlaps(
        #[case] region_chrom: &str,
        #[case] region_range: Option<(i32, i32)>,
        #[case] seqvar_chrom: &str,
        #[case] seqvar_start: i32,
        #[case] seqvar_end: i32,
        #[case] expected: bool,
    ) {
        let region = super::GenomicRegion {
            chrom: String::from(region_chrom),
            range: region_range.map(|(region_start, region_end)| super::Range {
                start: region_start,
                end: region_end,
            }),
        };
        assert_eq!(
            super::overlaps(&region, seqvar_chrom, seqvar_start, seqvar_end),
            expected
        );
    }
}
