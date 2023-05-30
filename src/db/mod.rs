//! Code supporting the `db *` sub commands.

use std::time::Instant;

pub mod compile;
pub mod conf;
pub mod genes;
pub mod mk_inhouse;
pub mod to_bin;

/// Compact column families of RocksDB instance.
pub fn rocksdb_compact_cf(
    cf_names: &[&str],
    db: &rocksdb::DBWithThreadMode<rocksdb::MultiThreaded>,
) -> Result<(), anyhow::Error> {
    let cfs_write = cf_names
        .iter()
        .map(|cf| db.cf_handle(cf).unwrap())
        .collect::<Vec<_>>();

    let mut compact_opt = rocksdb::CompactOptions::default();
    compact_opt.set_bottommost_level_compaction(rocksdb::BottommostLevelCompaction::Force);

    cfs_write
        .iter()
        .for_each(|cf| db.compact_range_cf_opt(cf, None::<&[u8]>, None::<&[u8]>, &compact_opt));
    let compaction_start = Instant::now();
    let mut last_printed = compaction_start;
    while db
        .property_int_value(rocksdb::properties::COMPACTION_PENDING)?
        .unwrap()
        > 0
        || db
            .property_int_value(rocksdb::properties::NUM_RUNNING_COMPACTIONS)?
            .unwrap()
            > 0
    {
        std::thread::sleep(std::time::Duration::from_millis(100));
        if last_printed.elapsed() > std::time::Duration::from_millis(1000) {
            log::info!(
                "... waiting for compaction for {:?}",
                compaction_start.elapsed()
            );
            last_printed = Instant::now();
        }
    }
    Ok(())
}

/// Tunin RocksDB options.
pub fn rocksdb_tuning(options: rocksdb::Options, wal_dir: Option<String>) -> rocksdb::Options {
    let mut options = options;

    options.create_if_missing(true);
    options.create_missing_column_families(true);

    options.prepare_for_bulk_load();

    // compress all files with Zstandard
    options.set_compression_per_level(&[]);
    options.set_compression_type(rocksdb::DBCompressionType::Zstd);
    // We only want to set level to 2 but have to set the rest as well using the
    // Rust interface. The (default) values for the other levels were taken from
    // the output of a RocksDB output folder created with default settings.
    options.set_compression_options(-14, 2, 0, 0);
    // options.set_zstd_max_train_bytes(100 * 1024);

    options.set_max_background_jobs(16);
    options.set_max_subcompactions(8);
    options.increase_parallelism(8);
    options.optimize_level_style_compaction(1 << 30);
    options.set_min_write_buffer_number(1);
    options.set_min_write_buffer_number_to_merge(1);
    options.set_write_buffer_size(1 << 30);
    options.set_target_file_size_base(1 << 30);
    options.set_compaction_style(rocksdb::DBCompactionStyle::Universal);

    if let Some(wal_dir) = wal_dir {
        options.set_wal_dir(wal_dir);
    }

    options.set_bottommost_compression_options(-14, 2, 0, 0, true);
    options.set_bottommost_compression_type(rocksdb::DBCompressionType::Zstd);

    options.set_disable_auto_compactions(true);
    options.optimize_for_point_lookup(1 << 26);

    options
}
