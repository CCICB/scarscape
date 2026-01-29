// use scarscape::io::manifest::ManifestReader;
//
// #[test]
// fn manifest_from_reader_parses_records_into_domain_structs() {
//     let data = b"sample_id\tsnv_vcf\tsv_vcf\tcnv_segments\n\
//                  S1\t/a/s1.snv.vcf.gz\t\t/a/s1.cnv.tsv\n\
//                  S2\t/a/s2.snv.vcf.gz\t/a/s2.sv.vcf.gz\t\n";
//
//     let it = ManifestReader::from_reader(&data[..]).unwrap();
//     let records: Vec<ManifestEntry> = it.map(|r| r.unwrap()).collect();
//
//     assert_eq!(records.len(), 2);
//     assert_eq!(records[0].sample_id.as_str(), "S1");
//     assert_eq!(records[0].snv_vcf.to_string_lossy(), "/a/s1.snv.vcf.gz");
//     assert!(records[0].sv_vcf.is_none());
//     assert_eq!(
//         records[0].cnv_segments.as_ref().unwrap().to_string_lossy(),
//         "/a/s1.cnv.tsv"
//     );
//
//     assert_eq!(records[1].sample_id.as_str(), "S2");
//     assert_eq!(
//         records[1].sv_vcf.as_ref().unwrap().to_string_lossy(),
//         "/a/s2.sv.vcf.gz"
//     );
//     assert!(records[1].cnv_segments.is_none());
// }
//
// #[test]
// fn manifest_missing_required_column_errors_immediately() {
//     let data = b"sample_id\tsv_vcf\nS1\t/a/s1.sv.vcf.gz\n";
//     let err = ManifestReader::from_reader(&data[..]).unwrap_err();
//     let msg = err.to_string();
//     assert!(msg.contains("missing required column"));
//     assert!(msg.contains("snv_vcf"));
// }
// // static manifest = include_str!("")
