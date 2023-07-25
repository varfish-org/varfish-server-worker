# Changelog

## [0.10.0](https://www.github.com/bihealth/varfish-server-worker/compare/v0.9.0...v0.10.0) (2023-07-25)


### Features

* write out HGNC ids for sv query results ([#150](https://www.github.com/bihealth/varfish-server-worker/issues/150)) ([#151](https://www.github.com/bihealth/varfish-server-worker/issues/151)) ([9e1074b](https://www.github.com/bihealth/varfish-server-worker/commit/9e1074b53572cf73344b3ac40e412d3d12f9f882))

## [0.9.0](https://www.github.com/bihealth/varfish-server-worker/compare/v0.8.0...v0.9.0) (2023-07-05)


### Features

* adding allelic balance for SVs ([#127](https://www.github.com/bihealth/varfish-server-worker/issues/127)) ([564a7ec](https://www.github.com/bihealth/varfish-server-worker/commit/564a7ec44754b461255a0a316ecd5ae59c66b002))
* allow filtering SVs based on effect on transcript ([#129](https://www.github.com/bihealth/varfish-server-worker/issues/129)) ([#134](https://www.github.com/bihealth/varfish-server-worker/issues/134)) ([1aab8c8](https://www.github.com/bihealth/varfish-server-worker/commit/1aab8c84af684c844c30dfd2a1cec72f0a3a2fc8))


### Bug Fixes

* dependabot config and updates ([#133](https://www.github.com/bihealth/varfish-server-worker/issues/133)) ([006e306](https://www.github.com/bihealth/varfish-server-worker/commit/006e3068500c3eb311cfa528e10bfbaf86046cb8))
* removing unused dependencies ([#124](https://www.github.com/bihealth/varfish-server-worker/issues/124)) ([2cee74b](https://www.github.com/bihealth/varfish-server-worker/commit/2cee74be09f97400ae3e59bc9a2bb435ff3956b2))
* specification of genomic range in query ([#126](https://www.github.com/bihealth/varfish-server-worker/issues/126)) ([6eddaef](https://www.github.com/bihealth/varfish-server-worker/commit/6eddaefc0b401d19755bd053ea79d09e31f4302a))

## [0.8.0](https://www.github.com/bihealth/varfish-server-worker/compare/v0.7.0...v0.8.0) (2023-06-28)


### ⚠ BREAKING CHANGES

* remove REST server related code (#122) (#123)
* replace db-compile with per-output commands (#112) (#121)

### Code Refactoring

* remove code that was moved to viguno/annonars ([#116](https://www.github.com/bihealth/varfish-server-worker/issues/116)) ([#117](https://www.github.com/bihealth/varfish-server-worker/issues/117)) ([50099dd](https://www.github.com/bihealth/varfish-server-worker/commit/50099dd01552f2b270728537ed34c01172de0e31))
* remove REST server related code ([#122](https://www.github.com/bihealth/varfish-server-worker/issues/122)) ([#123](https://www.github.com/bihealth/varfish-server-worker/issues/123)) ([0b68809](https://www.github.com/bihealth/varfish-server-worker/commit/0b68809275294dd295b72fc73faa8fc26c095723))
* replace db-compile with per-output commands ([#112](https://www.github.com/bihealth/varfish-server-worker/issues/112)) ([#121](https://www.github.com/bihealth/varfish-server-worker/issues/121)) ([5da2799](https://www.github.com/bihealth/varfish-server-worker/commit/5da27991765cdc294a4114aa57dbb2464e112d6e))

## [0.7.0](https://www.github.com/bihealth/varfish-server-worker/compare/v0.6.1...v0.7.0) (2023-06-09)


### Features

* add "db genes query" command ([#109](https://www.github.com/bihealth/varfish-server-worker/issues/109)) ([27e995c](https://www.github.com/bihealth/varfish-server-worker/commit/27e995c9a1127f2e2fc41a5481b2ae07a3518ea7))
* adding "pheno prune" command against all-parents annotation ([#100](https://www.github.com/bihealth/varfish-server-worker/issues/100)) ([#114](https://www.github.com/bihealth/varfish-server-worker/issues/114)) ([c77b2d4](https://www.github.com/bihealth/varfish-server-worker/commit/c77b2d4dbd2a3708368fccc1371b2f6364042440))
* adding genes database ([#85](https://www.github.com/bihealth/varfish-server-worker/issues/85)) ([#90](https://www.github.com/bihealth/varfish-server-worker/issues/90)) ([634e281](https://www.github.com/bihealth/varfish-server-worker/commit/634e281c2528f3a6b072e36618509faa74c796e7))
* adding sequence variant databases ([#94](https://www.github.com/bihealth/varfish-server-worker/issues/94)) ([74f985a](https://www.github.com/bihealth/varfish-server-worker/commit/74f985ae6ba945d2cda1c4d5ceab277d309ba675))
* individual paths for "db genes build" ([#106](https://www.github.com/bihealth/varfish-server-worker/issues/106)) ([0a10a4e](https://www.github.com/bihealth/varfish-server-worker/commit/0a10a4e0f93129fe01cb9df77b80272a9652dc23))
* move from flatbuffers to protocol buffers ([#89](https://www.github.com/bihealth/varfish-server-worker/issues/89)) ([#99](https://www.github.com/bihealth/varfish-server-worker/issues/99)) ([cff7ace](https://www.github.com/bihealth/varfish-server-worker/commit/cff7ace317d705fe46fa106ea3a56a1fe9b34457))
* support for dbNSFP gene information ([#97](https://www.github.com/bihealth/varfish-server-worker/issues/97)) ([#102](https://www.github.com/bihealth/varfish-server-worker/issues/102)) ([0448bcb](https://www.github.com/bihealth/varfish-server-worker/commit/0448bcba3e39eca1d69c3553d711feba9e275b46))
* use protobuf/prost for gene information storage ([#96](https://www.github.com/bihealth/varfish-server-worker/issues/96)) ([#101](https://www.github.com/bihealth/varfish-server-worker/issues/101)) ([6a97fb6](https://www.github.com/bihealth/varfish-server-worker/commit/6a97fb6c58c9800c519969c99b0373b915420e72))

### [0.6.1](https://www.github.com/bihealth/varfish-server-worker/compare/v0.6.0...v0.6.1) (2023-05-03)


### Bug Fixes

* "db mk-inhouse" command for breakends ([#79](https://www.github.com/bihealth/varfish-server-worker/issues/79)) ([9b33e1d](https://www.github.com/bihealth/varfish-server-worker/commit/9b33e1d31cdbef31c636da40d1394eac63fb322f))
* command "sv query" for types BND and INS ([#81](https://www.github.com/bihealth/varfish-server-worker/issues/81)) ([5be2274](https://www.github.com/bihealth/varfish-server-worker/commit/5be2274c86f86f3ae32b255f3e778ddebc1c7a28))
* correctly interpret "callers" column ([#82](https://www.github.com/bihealth/varfish-server-worker/issues/82)) ([8e083c8](https://www.github.com/bihealth/varfish-server-worker/commit/8e083c87b3d5e7d263ade2d4540f06250b20bc6b))
* overlap queries for BND and INS ([#83](https://www.github.com/bihealth/varfish-server-worker/issues/83)) ([4de4eda](https://www.github.com/bihealth/varfish-server-worker/commit/4de4eda55c2caa34881305339a194eae2125b2f1))

## [0.6.0](https://www.github.com/bihealth/varfish-server-worker/compare/v0.5.1...v0.6.0) (2023-05-02)


### Features

* add "serve pheno" REST API service ([#72](https://www.github.com/bihealth/varfish-server-worker/issues/72)) ([#73](https://www.github.com/bihealth/varfish-server-worker/issues/73)) ([17c9dd5](https://www.github.com/bihealth/varfish-server-worker/commit/17c9dd53c0a14949087e03b840dd0e7f50a1b55e))

### [0.5.1](https://www.github.com/bihealth/varfish-server-worker/compare/v0.5.0...v0.5.1) (2023-03-14)


### Bug Fixes

* reciprocal overlap computation ([#69](https://www.github.com/bihealth/varfish-server-worker/issues/69)) ([6165d9f](https://www.github.com/bihealth/varfish-server-worker/commit/6165d9f4f8df412ccb5bfbd6a8d4783f6f3d0e15))
* region deserialization in sv query ([#70](https://www.github.com/bihealth/varfish-server-worker/issues/70)) ([7a577e3](https://www.github.com/bihealth/varfish-server-worker/commit/7a577e3f5f2711f0247e6fc0e7f332af450e1197))

## [0.5.0](https://www.github.com/bihealth/varfish-server-worker/compare/v0.4.0...v0.5.0) (2023-02-13)


### Features

* annotate SVs with masked sequence overlap ([#63](https://www.github.com/bihealth/varfish-server-worker/issues/63)) ([#66](https://www.github.com/bihealth/varfish-server-worker/issues/66)) ([0fd9403](https://www.github.com/bihealth/varfish-server-worker/commit/0fd940362f820913a71fbca5c53194206e097792))
* implement gene allow lists for "sv query" ([#65](https://www.github.com/bihealth/varfish-server-worker/issues/65)) ([#67](https://www.github.com/bihealth/varfish-server-worker/issues/67)) ([aa9da21](https://www.github.com/bihealth/varfish-server-worker/commit/aa9da21b27a0ba30b8a2589f543a44643815abc2))
* serve genome browser tracks from database ([#59](https://www.github.com/bihealth/varfish-server-worker/issues/59)) ([#60](https://www.github.com/bihealth/varfish-server-worker/issues/60)) ([040c522](https://www.github.com/bihealth/varfish-server-worker/commit/040c5227b251d76594acfa93f714305e0910eca5))
* write out effective and matched genotype from rules ([#64](https://www.github.com/bihealth/varfish-server-worker/issues/64)) ([#68](https://www.github.com/bihealth/varfish-server-worker/issues/68)) ([2c2c793](https://www.github.com/bihealth/varfish-server-worker/commit/2c2c7934d3fe333d72b709847a1c4c4a48375973))


### Bug Fixes

* hgnc xlink path ([#62](https://www.github.com/bihealth/varfish-server-worker/issues/62)) ([985cc33](https://www.github.com/bihealth/varfish-server-worker/commit/985cc334e99e6d693c3a665b34a825d28f4c81b6))

## [0.4.0](https://www.github.com/bihealth/varfish-server-worker/compare/v0.3.0...v0.4.0) (2023-02-09)


### Features

* add "server rest" sub command ([#57](https://www.github.com/bihealth/varfish-server-worker/issues/57)) ([#58](https://www.github.com/bihealth/varfish-server-worker/issues/58)) ([4b41010](https://www.github.com/bihealth/varfish-server-worker/commit/4b41010fb4acddf18f44904358d61a21dcfe8df6))
* implement "db compile" sub command ([#54](https://www.github.com/bihealth/varfish-server-worker/issues/54)) ([#55](https://www.github.com/bihealth/varfish-server-worker/issues/55)) ([7aa8c69](https://www.github.com/bihealth/varfish-server-worker/commit/7aa8c692bcef4d39c7b8aa6b82be5a5087ea8658))

## [0.3.0](https://www.github.com/bihealth/varfish-server-worker/compare/v0.2.1...v0.3.0) (2023-02-01)


### ⚠ BREAKING CHANGES

* enums serialization case and query issues (#41)

### Features

* allow configuring whether missing call_info entry is pass ([#47](https://www.github.com/bihealth/varfish-server-worker/issues/47)) ([#49](https://www.github.com/bihealth/varfish-server-worker/issues/49)) ([4d0f64d](https://www.github.com/bihealth/varfish-server-worker/commit/4d0f64da21eb883cd7c31d56c9eacddd9fb4361b))
* make clinvar/pathogenic SV overlap configurable ([#39](https://www.github.com/bihealth/varfish-server-worker/issues/39)) ([af36bd4](https://www.github.com/bihealth/varfish-server-worker/commit/af36bd47dc50740543b6125e555aad7fe25e81e4))
* make TAD set configurable in query ([#51](https://www.github.com/bihealth/varfish-server-worker/issues/51)) ([82d4995](https://www.github.com/bihealth/varfish-server-worker/commit/82d4995ca727e46018a30520da48b31d774b7e23))
* removing minimal overlap for pathogenic variants again ([#53](https://www.github.com/bihealth/varfish-server-worker/issues/53)) ([3a5d427](https://www.github.com/bihealth/varfish-server-worker/commit/3a5d427ac3e6feaac74ee555752b56faf04ffd27))
* write out bg db overlap counts in "sv query" ([#44](https://www.github.com/bihealth/varfish-server-worker/issues/44)) ([#46](https://www.github.com/bihealth/varfish-server-worker/issues/46)) ([7286d0d](https://www.github.com/bihealth/varfish-server-worker/commit/7286d0d52d51eae39145772f9121f4c960f3fdd2))
* write out distance to next TAD boundary ([#45](https://www.github.com/bihealth/varfish-server-worker/issues/45)) ([#50](https://www.github.com/bihealth/varfish-server-worker/issues/50)) ([5ef72f3](https://www.github.com/bihealth/varfish-server-worker/commit/5ef72f32436c3e7d2116a26193935c7ed0db7b5a))
* writing out payload.sv_length ([#43](https://www.github.com/bihealth/varfish-server-worker/issues/43)) ([368621d](https://www.github.com/bihealth/varfish-server-worker/commit/368621dfbebfb1e0b1c1eb40f466277780cca82e))


### Bug Fixes

* enums serialization case and query issues ([#41](https://www.github.com/bihealth/varfish-server-worker/issues/41)) ([7cd73aa](https://www.github.com/bihealth/varfish-server-worker/commit/7cd73aae69421998500571decd9417f026824135))


### Miscellaneous Chores

* set next 0.3.0 ([101615d](https://www.github.com/bihealth/varfish-server-worker/commit/101615dbdaf3c302b226aaaec77866ba31aca5ae))

### [0.2.1](https://www.github.com/bihealth/varfish-server-worker/compare/v0.2.0...v0.2.1) (2023-01-31)


### Bug Fixes

* pinning flatc-rust ([#37](https://www.github.com/bihealth/varfish-server-worker/issues/37)) ([e9d1d89](https://www.github.com/bihealth/varfish-server-worker/commit/e9d1d89a16fe5c962cfe97e215f378ba138968e4))

## [0.2.0](https://www.github.com/bihealth/varfish-server-worker/compare/v0.1.0...v0.2.0) (2023-01-31)


### Features

* adding --max-results for "sv query" sub command ([#32](https://www.github.com/bihealth/varfish-server-worker/issues/32)) ([#33](https://www.github.com/bihealth/varfish-server-worker/issues/33)) ([4aa7094](https://www.github.com/bihealth/varfish-server-worker/commit/4aa709468d0e8a3bdd3aaa78b6376b9a31f08c25))
* annotate ACMG and OMIM genes in "sv query" ([#31](https://www.github.com/bihealth/varfish-server-worker/issues/31)) ([#35](https://www.github.com/bihealth/varfish-server-worker/issues/35)) ([9058e2a](https://www.github.com/bihealth/varfish-server-worker/commit/9058e2aca09e97e908009af7f00c3e6663a73ee6))

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
