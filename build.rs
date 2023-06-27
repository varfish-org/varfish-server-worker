// The custom build script, needed as we use protocolbuffers.

fn main() {
    println!("cargo:rerun-if-changed=varfish/v1/clinvar.proto");
    println!("cargo:rerun-if-changed=varfish/v1/sv.proto");
    prost_build::Config::new()
        .protoc_arg("-Isrc/proto")
        // Add serde serialization and deserialization to the generated code.
        .type_attribute(".", "#[derive(serde::Serialize, serde::Deserialize)]")
        // Skip serializing `None` values.
        .type_attribute(".", "#[serde_with::skip_serializing_none]")
        // Define the protobuf files to compile.
        .compile_protos(
            &["varfish/v1/clinvar.proto", "varfish/v1/sv.proto"],
            &["src/"],
        )
        .unwrap();
}
