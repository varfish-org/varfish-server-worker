//! Definition of records used by the ``sv`` sub command.

use std::{fs, io, path::Path};

use anyhow::anyhow;
use md5::{Digest, Md5};
use serde::{Deserialize, Serialize};
use sha2::Sha256;
use std::io::Read;
use tracing::debug;

/// Configuration for the database backing the SV annotation.
#[derive(Serialize, Deserialize, PartialEq, Debug)]
pub struct DbConf {
    /// Configuration of public databases.
    pub public_dbs: PublicDbsConf,
    /// Configuration of TAD boundaries.
    pub tads: TadsConf,
    /// Configuration of regulatory features.
    pub regulatory: RegulatoryConf,
    /// Gene configuration
    pub genes: GenesConf,
}

/// Configuration for public databases.
#[derive(Serialize, Deserialize, PartialEq, Debug)]
pub struct PublicDbsConf {
    /// Radius around BND sites used when building the database.
    #[serde(default = "default_slack_bnd")]
    pub slack_bnd: u32,
    /// Radius around INS sites used when building the database.
    #[serde(default = "default_slack_ins")]
    pub slack_ins: u32,
    /// Minimal reciprocal overlap for SVs of the same type, used when building the database.
    #[serde(default = "default_min_overlap")]
    pub min_overlap: f32,

    /// Relative path to the file with gnomAD-SV database with checksum.
    pub gnomad_sv: PathAndChecksum,
    /// Relative path to the file with dbVar SV database with checksum.
    pub dbvar: PathAndChecksum,
    /// Relative path to the file with DGV database with checksum.
    pub dgv: PathAndChecksum,
    /// Relative path to the file with DGV GS database with checksum.
    pub dgv_gs: PathAndChecksum,
    /// Relative path to the file with ExAC database with checksum.
    pub exac: PathAndChecksum,
}

fn default_slack_bnd() -> u32 {
    50
}

fn default_slack_ins() -> u32 {
    50
}

fn default_min_overlap() -> f32 {
    0.8
}

/// A relative path to a file with its checksum.
#[derive(Serialize, Deserialize, PartialEq, Debug)]
pub struct PathAndChecksum {
    /// The (relative) path.
    pub path: String,
    /// Optional MD5 checksum.
    pub md5: Option<String>,
    /// Optional SHA256 checksum.
    pub sha256: Option<String>,
}

/// TAD configuration.
#[derive(Serialize, Deserialize, PartialEq, Debug)]
pub struct TadsConf {
    /// Maximal distance to TAD to display for.
    #[serde(default = "default_max_dist")]
    pub max_dist: u32,
    /// Path and checksum for Dixon hESC
    pub hesc: PathAndChecksum,
    /// Path and checksum for Dixon IMR90
    pub imr90: PathAndChecksum,
}

fn default_max_dist() -> u32 {
    10_000
}

/// Configuration of regulatory features.
#[derive(Serialize, Deserialize, PartialEq, Debug)]
pub struct RegulatoryConf {
    /// Path and checksum for ENSEMBL regulatory build
    pub ensembl: PathAndChecksum,
    /// Path and checksum for VISTA validations
    pub vista: PathAndChecksum,
}

/// Configuration of genes.
#[derive(Serialize, Deserialize, PartialEq, Debug)]
pub struct GenesConf {
    /// Gene identifier interlink table.
    pub xlink: PathAndChecksum,
    /// ACMG incidental finding list information
    pub acmg: PathAndChecksum,
    /// Annotation of genes with associated HPO terms and OMIM diseases from NCBI medgen.
    pub medgen: PathAndChecksum,

    /// RefSeq gene configuration
    pub refseq: GenesDetailConf,
    /// ENSEMBL gene configuration
    pub ensembl: GenesDetailConf,
}

/// Configuration of gene information database.
#[derive(Serialize, Deserialize, PartialEq, Debug)]
pub struct GenesDetailConf {
    /// Regions for overlapping
    pub regions: PathAndChecksum,
    /// RocksDB table indexed by gene ID with annotation payload JSON
    pub anno_json: PathAndChecksum,
}

/// Check whether the given path is valid with the given path to the configuration.
///
/// Return `Vec` with error messages or `None` if everything validates.  In case of any errors
/// (validation failure is a result, not an error), an `anyhow::Error` is returned.
pub fn sanity_check_db(
    path_db: &Path,
    path_config: &Path,
    checksums: bool,
) -> Result<Option<Vec<String>>, anyhow::Error> {
    let mut result = Vec::new();

    let toml_str = match std::fs::read_to_string(path_config) {
        Err(err) => {
            return Err(anyhow!(
                "Could not open config file {:?}: {}",
                path_config,
                err
            ))
        }
        Ok(toml_str) => toml_str,
    };
    let conf: DbConf = match toml::from_str(&toml_str) {
        Err(err) => {
            return Err(anyhow!(
                "Could not load configuration file {:?}: {}",
                path_config,
                err
            ))
        }
        Ok(conf) => conf,
    };

    check_path_and_checksum(path_db, &conf.public_dbs.gnomad_sv, checksums, &mut result)?;
    check_path_and_checksum(path_db, &conf.public_dbs.dbvar, checksums, &mut result)?;
    check_path_and_checksum(path_db, &conf.public_dbs.dgv, checksums, &mut result)?;
    check_path_and_checksum(path_db, &conf.public_dbs.dgv_gs, checksums, &mut result)?;
    check_path_and_checksum(path_db, &conf.public_dbs.exac, checksums, &mut result)?;
    check_path_and_checksum(path_db, &conf.tads.hesc, checksums, &mut result)?;
    check_path_and_checksum(path_db, &conf.tads.imr90, checksums, &mut result)?;
    check_path_and_checksum(path_db, &conf.tads.imr90, checksums, &mut result)?;
    check_path_and_checksum(path_db, &conf.regulatory.ensembl, checksums, &mut result)?;
    check_path_and_checksum(path_db, &conf.regulatory.vista, checksums, &mut result)?;
    check_path_and_checksum(path_db, &conf.genes.xlink, checksums, &mut result)?;
    check_path_and_checksum(path_db, &conf.genes.acmg, checksums, &mut result)?;
    check_path_and_checksum(path_db, &conf.genes.refseq.regions, checksums, &mut result)?;
    check_path_and_checksum(
        path_db,
        &conf.genes.refseq.anno_json,
        checksums,
        &mut result,
    )?;
    check_path_and_checksum(path_db, &conf.genes.ensembl.regions, checksums, &mut result)?;
    check_path_and_checksum(
        path_db,
        &conf.genes.ensembl.anno_json,
        checksums,
        &mut result,
    )?;

    if result.is_empty() {
        Ok(None)
    } else {
        Ok(Some(result))
    }
}

/// Check that the `path_and_checksum` relative to `base_path` exist and the checksum
/// is correct (if any).
fn check_path_and_checksum(
    base_path: &Path,
    path_and_checksum: &PathAndChecksum,
    checksums: bool,
    result: &mut Vec<String>,
) -> Result<(), anyhow::Error> {
    let joined_path = base_path.join(&path_and_checksum.path);
    let full_path = if path_and_checksum.path.starts_with('/') {
        Path::new(&path_and_checksum.path)
    } else {
        joined_path.as_path()
    };

    if !full_path.exists() {
        result.push(format!("file {:?} does not exist", &full_path));
        return Ok(());
    }

    if checksums {
        if let Some(sha256) = &path_and_checksum.sha256 {
            debug!(
                "SHA256 checksum verification for {}",
                &path_and_checksum.path
            );
            let mut file = fs::File::open(full_path)?;
            let mut hasher = Sha256::new();
            let n = io::copy(&mut file, &mut hasher)?;
            let hash = hasher.finalize();
            let mut buf = [0u8; 64];
            let checksum = base16ct::lower::encode_str(&hash, &mut buf).unwrap();
            if checksum != sha256 {
                result.push(format!(
                    "file {:?} (with {} bytes) has checksum sha256:{} instead of sha256:{}",
                    &full_path, n, &checksum, &sha256
                ));
            }
        } else {
            debug!("SHA256 checksum unknown for {}", &path_and_checksum.path);
        }

        if let Some(md5) = &path_and_checksum.md5 {
            debug!("MD5 checksum verification for {}", &path_and_checksum.path);
            let mut file = fs::File::open(full_path)?;
            let mut hasher = Md5::new();
            const BUF_LEN: usize = 65_536;
            let mut bytes_read = 0;
            let mut buffer = [0; BUF_LEN];

            loop {
                let n = file.read(&mut buffer)?;
                hasher.update(&buffer[..n]);
                bytes_read += n;
                if n != BUF_LEN {
                    break;
                }
            }
            let hash = hasher.finalize();
            let mut buf = [0u8; 64];
            let checksum = base16ct::lower::encode_str(&hash, &mut buf).unwrap();
            if checksum != md5 {
                result.push(format!(
                    "file {:?} (with {} bytes) has checksum md5:{} instead of md5:{}",
                    &full_path, bytes_read, &checksum, &md5
                ));
            }
        } else {
            debug!("MD5> checksum unknown for {}", &path_and_checksum.path);
        }
    } else {
        debug!(
            "checksums verification disabled for {}",
            &path_and_checksum.path
        );
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use pretty_assertions::assert_eq;

    use crate::sv::conf::{
        sanity_check_db, GenesConf, GenesDetailConf, PathAndChecksum, PublicDbsConf,
        RegulatoryConf, TadsConf,
    };

    use super::DbConf;

    #[test]
    fn test_check_path_and_checksum_good() -> Result<(), anyhow::Error> {
        let errors = sanity_check_db(
            &Path::new("tests/sv/conf/good"),
            &Path::new("tests/sv/conf/full.toml"),
            true,
        )?;
        assert!(errors.is_none());

        Ok(())
    }

    #[test]
    fn test_check_path_and_checksum_bad() -> Result<(), anyhow::Error> {
        let errors = sanity_check_db(
            &Path::new("tests/sv/conf/bad"),
            &Path::new("tests/sv/conf/full.toml"),
            true,
        )?;
        assert!(errors.is_some());
        assert_eq!(errors.unwrap().len(), 32);

        Ok(())
    }

    #[test]
    fn test_check_path_and_checksum_nonexisting() -> Result<(), anyhow::Error> {
        let errors = sanity_check_db(
            &Path::new("tests/sv/conf/nonexisting"),
            &Path::new("tests/sv/conf/full.toml"),
            true,
        )?;
        assert!(errors.is_some());
        assert_eq!(errors.unwrap().len(), 16);

        Ok(())
    }

    #[test]
    fn test_parse_config_full() -> Result<(), anyhow::Error> {
        let toml_path = "tests/sv/conf/full.toml";
        let toml_str = std::fs::read_to_string(toml_path)?;
        let toml_data: DbConf = toml::from_str(&toml_str)?;

        assert_eq!(
            toml_data,
            DbConf {
                public_dbs: PublicDbsConf {
                    slack_bnd: 50,
                    slack_ins: 50,
                    min_overlap: 0.8,
                    gnomad_sv: PathAndChecksum {
                        path: "public/gnomad_sv.bin".to_owned(),
                        md5: Some("d41d8cd98f00b204e9800998ecf8427e".to_owned()),
                        sha256: Some(
                            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
                                .to_owned()
                        ),
                    },
                    dbvar: PathAndChecksum {
                        path: "public/dbvar.bin".to_owned(),
                        md5: Some("d41d8cd98f00b204e9800998ecf8427e".to_owned()),
                        sha256: Some(
                            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
                                .to_owned()
                        )
                    },
                    dgv: PathAndChecksum {
                        path: "public/dgv.bin".to_owned(),
                        md5: Some("d41d8cd98f00b204e9800998ecf8427e".to_owned()),
                        sha256: Some(
                            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
                                .to_owned()
                        ),
                    },
                    dgv_gs: PathAndChecksum {
                        path: "public/dgv_gs.bin".to_owned(),
                        md5: Some("d41d8cd98f00b204e9800998ecf8427e".to_owned()),
                        sha256: Some(
                            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
                                .to_owned()
                        ),
                    },
                    exac: PathAndChecksum {
                        path: "public/exac.bin".to_owned(),
                        md5: Some("d41d8cd98f00b204e9800998ecf8427e".to_owned()),
                        sha256: Some(
                            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
                                .to_owned()
                        ),
                    }
                },
                tads: TadsConf {
                    max_dist: 10_000,
                    hesc: PathAndChecksum {
                        path: "tads/hesc.bin".to_owned(),
                        md5: Some("d41d8cd98f00b204e9800998ecf8427e".to_owned()),
                        sha256: Some(
                            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
                                .to_owned()
                        ),
                    },
                    imr90: PathAndChecksum {
                        path: "tads/imr90.bin".to_owned(),
                        md5: Some("d41d8cd98f00b204e9800998ecf8427e".to_owned()),
                        sha256: Some(
                            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
                                .to_owned()
                        ),
                    }
                },
                regulatory: RegulatoryConf {
                    ensembl: PathAndChecksum {
                        path: "regulatory/ensembl.bin".to_owned(),
                        md5: Some("d41d8cd98f00b204e9800998ecf8427e".to_owned()),
                        sha256: Some(
                            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
                                .to_owned()
                        ),
                    },
                    vista: PathAndChecksum {
                        path: "regulatory/vista.bin".to_owned(),
                        md5: Some("d41d8cd98f00b204e9800998ecf8427e".to_owned()),
                        sha256: Some(
                            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
                                .to_owned()
                        ),
                    },
                },
                genes: GenesConf {
                    xlink: PathAndChecksum {
                        path: "genes/xlink.bin".to_owned(),
                        md5: Some("d41d8cd98f00b204e9800998ecf8427e".to_owned()),
                        sha256: Some(
                            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
                                .to_owned()
                        ),
                    },
                    acmg: PathAndChecksum {
                        path: "genes/acmg.json".to_owned(),
                        md5: Some("d41d8cd98f00b204e9800998ecf8427e".to_owned()),
                        sha256: Some(
                            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
                                .to_owned()
                        ),
                    },
                    medgen: PathAndChecksum {
                        path: "genes/medgen.json".to_owned(),
                        md5: Some("d41d8cd98f00b204e9800998ecf8427e".to_owned()),
                        sha256: Some(
                            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
                                .to_owned()
                        ),
                    },
                    refseq: GenesDetailConf {
                        regions: PathAndChecksum {
                            path: "genes/refseq/regions.bin".to_owned(),
                            md5: Some("d41d8cd98f00b204e9800998ecf8427e".to_owned()),
                            sha256: Some(
                                "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
                                    .to_owned()
                            ),
                        },
                        anno_json: PathAndChecksum {
                            path: "genes/refseq/anno_json.db".to_owned(),
                            md5: Some("d41d8cd98f00b204e9800998ecf8427e".to_owned()),
                            sha256: Some(
                                "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
                                    .to_owned()
                            ),
                        }
                    },
                    ensembl: GenesDetailConf {
                        regions: PathAndChecksum {
                            path: "genes/ensembl/regions.bin".to_owned(),
                            md5: Some("d41d8cd98f00b204e9800998ecf8427e".to_owned()),
                            sha256: Some(
                                "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
                                    .to_owned()
                            ),
                        },
                        anno_json: PathAndChecksum {
                            path: "genes/ensembl/anno_json.db".to_owned(),
                            md5: Some("d41d8cd98f00b204e9800998ecf8427e".to_owned()),
                            sha256: Some(
                                "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
                                    .to_owned()
                            ),
                        }
                    },
                }
            }
        );

        Ok(())
    }

    #[test]
    fn test_parse_config_minimal() -> Result<(), anyhow::Error> {
        let toml_path = "tests/sv/conf/minimal.toml";
        let toml_str = std::fs::read_to_string(toml_path)?;
        let toml_data: DbConf = toml::from_str(&toml_str)?;

        assert_eq!(
            toml_data,
            DbConf {
                public_dbs: PublicDbsConf {
                    slack_bnd: 50,
                    slack_ins: 50,
                    min_overlap: 0.8,
                    gnomad_sv: PathAndChecksum {
                        path: "public/gnomad_sv.bin".to_owned(),
                        md5: None,
                        sha256: None,
                    },
                    dbvar: PathAndChecksum {
                        path: "public/dbvar.bin".to_owned(),
                        md5: None,
                        sha256: None,
                    },
                    dgv: PathAndChecksum {
                        path: "public/dgv.bin".to_owned(),
                        md5: None,
                        sha256: None,
                    },
                    dgv_gs: PathAndChecksum {
                        path: "public/dgv_gs.bin".to_owned(),
                        md5: None,
                        sha256: None,
                    },
                    exac: PathAndChecksum {
                        path: "public/exac.bin".to_owned(),
                        md5: None,
                        sha256: None,
                    }
                },
                tads: TadsConf {
                    max_dist: 10_000,
                    hesc: PathAndChecksum {
                        path: "tads/hesc.bin".to_owned(),
                        md5: None,
                        sha256: None,
                    },
                    imr90: PathAndChecksum {
                        path: "tads/imr90.bin".to_owned(),
                        md5: None,
                        sha256: None,
                    }
                },
                regulatory: RegulatoryConf {
                    ensembl: PathAndChecksum {
                        path: "regulatory/ensembl.bin".to_owned(),
                        md5: None,
                        sha256: None,
                    },
                    vista: PathAndChecksum {
                        path: "regulatory/vista.bin".to_owned(),
                        md5: None,
                        sha256: None,
                    },
                },
                genes: GenesConf {
                    xlink: PathAndChecksum {
                        path: "genes/xlink.bin".to_owned(),
                        md5: None,
                        sha256: None,
                    },
                    acmg: PathAndChecksum {
                        path: "genes/acmg.json".to_owned(),
                        md5: None,
                        sha256: None,
                    },
                    medgen: PathAndChecksum {
                        path: "genes/medgen.json".to_owned(),
                        md5: None,
                        sha256: None,
                    },
                    refseq: GenesDetailConf {
                        regions: PathAndChecksum {
                            path: "genes/refseq/regions.bin".to_owned(),
                            md5: None,
                            sha256: None,
                        },
                        anno_json: PathAndChecksum {
                            path: "genes/refseq/anno_json.db".to_owned(),
                            md5: None,
                            sha256: None,
                        }
                    },
                    ensembl: GenesDetailConf {
                        regions: PathAndChecksum {
                            path: "genes/ensembl/regions.bin".to_owned(),
                            md5: None,
                            sha256: None,
                        },
                        anno_json: PathAndChecksum {
                            path: "genes/ensembl/anno_json.db".to_owned(),
                            md5: None,
                            sha256: None,
                        }
                    },
                }
            }
        );

        Ok(())
    }
}
