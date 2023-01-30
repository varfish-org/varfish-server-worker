# Changelog

## 0.1.0 (2023-01-30)


### Features

* "sv query" now writes out TSV with information for a MVP ([#28](https://www.github.com/bihealth/varfish-server-worker/issues/28)) ([6a928c8](https://www.github.com/bihealth/varfish-server-worker/commit/6a928c8d92d136bc5f53d8204aae3da85815c07f))
* add "sv build-bgdb" command ([#18](https://www.github.com/bihealth/varfish-server-worker/issues/18)) ([#19](https://www.github.com/bihealth/varfish-server-worker/issues/19)) ([00cdf41](https://www.github.com/bihealth/varfish-server-worker/commit/00cdf413ddd9b18994a8ba4f0e4902aa37adbaa1))
* annotating SVs with ClinVar VCV ([#27](https://www.github.com/bihealth/varfish-server-worker/issues/27)) ([14af2bc](https://www.github.com/bihealth/varfish-server-worker/commit/14af2bc7c41126a43639fd75b4a277bf2c8ceddd))
* bootstrapping sv-query sub command ([#2](https://www.github.com/bihealth/varfish-server-worker/issues/2)) ([9cb0c95](https://www.github.com/bihealth/varfish-server-worker/commit/9cb0c95a390a99293bee6862ef80efc3b9437429))
* implement proper counting of background records ([#13](https://www.github.com/bihealth/varfish-server-worker/issues/13)) ([#14](https://www.github.com/bihealth/varfish-server-worker/issues/14)) ([a0965cc](https://www.github.com/bihealth/varfish-server-worker/commit/a0965cc1eeb5c73c6becc9f7c1d62f4d12030093))
* implement public/in-house SV database interval tree construction ([#4](https://www.github.com/bihealth/varfish-server-worker/issues/4)) ([#5](https://www.github.com/bihealth/varfish-server-worker/issues/5)) ([03442aa](https://www.github.com/bihealth/varfish-server-worker/commit/03442aa89769a12eb75fbad01d791bb9e29e5445))
* implement query schemas (10) ([#11](https://www.github.com/bihealth/varfish-server-worker/issues/11)) ([cc087d7](https://www.github.com/bihealth/varfish-server-worker/commit/cc087d7d8198d901feb97745f18627efcc391c3c))
* implementing PoC query feature with SV VCF ([#6](https://www.github.com/bihealth/varfish-server-worker/issues/6)) ([#7](https://www.github.com/bihealth/varfish-server-worker/issues/7)) ([925ebd8](https://www.github.com/bihealth/varfish-server-worker/commit/925ebd85bf8d50d9debd6a5dca9e8851bc90955a))
* make "sv query-next" feature-complete for a MVP ([#24](https://www.github.com/bihealth/varfish-server-worker/issues/24)) ([#25](https://www.github.com/bihealth/varfish-server-worker/issues/25)) ([b7c873c](https://www.github.com/bihealth/varfish-server-worker/commit/b7c873c3f46c8bb393a26407037fe9f90995ec63))
* PoC implementation of variant filtration ([#15](https://www.github.com/bihealth/varfish-server-worker/issues/15)) ([#16](https://www.github.com/bihealth/varfish-server-worker/issues/16)) ([12201ca](https://www.github.com/bihealth/varfish-server-worker/commit/12201ca6f3ad1db7564252ffa905fc438895d3c0))
* read query from JSON file ([#22](https://www.github.com/bihealth/varfish-server-worker/issues/22)) ([#23](https://www.github.com/bihealth/varfish-server-worker/issues/23)) ([4e41fbc](https://www.github.com/bihealth/varfish-server-worker/commit/4e41fbc72eff9f25d7cd59b714c89a369e204eb4))
* switching CLI to multiple sub command levels ([#17](https://www.github.com/bihealth/varfish-server-worker/issues/17)) ([9af47df](https://www.github.com/bihealth/varfish-server-worker/commit/9af47dfc358b70a9fcf7e5bcb50a1e2b14e4fa0c))
* using reasonble settings for sv query poc ([#20](https://www.github.com/bihealth/varfish-server-worker/issues/20)) ([#21](https://www.github.com/bihealth/varfish-server-worker/issues/21)) ([26c0ba0](https://www.github.com/bihealth/varfish-server-worker/commit/26c0ba021d21d0969cfc68806a136060902c273b))


### Performance Improvements

* adjust interval types ([#6](https://www.github.com/bihealth/varfish-server-worker/issues/6)) ([#9](https://www.github.com/bihealth/varfish-server-worker/issues/9)) ([76d254b](https://www.github.com/bihealth/varfish-server-worker/commit/76d254be54b499455a1e02d0bde6af08aea9913f))