# 植生被度階級データの状態空間モデリング

第132回日本森林学会大会の企画シンポジウム「階層モデリングは森林の隠れた生態的状態の推測を可能にする」における講演「植生被度階級データの状態空間モデリング」の発表資料および実行コード (R Markdown形式)

## 発表資料

- 植生被度階級データの状態空間モデリング（準備中）

## ファイル

- jfs132.Rmd: 発表資料のR Markdownファイル
    - cover.stan: 被度階級データにベータ分布をあてはめるStanモデルファイル
    - cover_ssm.stan: 時系列の被度階級データの状態空間モデルのStanモデルファイル
    - cover_ssm2.stan: 過程モデルをトレンドモデルにしたStanモデルファイル

## 参考文献

- Damgaard C. (2012) Trend analyses of hierarchical pin-point cover data. Ecology 93:1269–1274. doi:[10.1890/11-1499.1](https://doi.org/10.1890/11-1499.1)
- Damgaard C.F., Irvine K.M. (2019) Using the beta distribution to analyse plant cover data. Journal of Ecology 107:2747–2759. [doi:10.1111/1365-2745.13200](https://doi.org/10.1111/1365-2745.13200)
- Herpigny B., Gosselin F. (2015) Analyzing plant cover class data quantitatively: Customized zero-inflated cumulative beta distributions show promising results. Ecological Informatics 26:18–26. [doi:10.1016/j.ecoinf.2014.12.002](https://doi.org/10.1016/j.ecoinf.2014.12.002)
- Irvine K.M., Rodhouse T.J., Keren I.N. (2016) Extending ordinal regression with a latent zero-augmented beta distribution. Journal of Agricultural, Biological and Environmental Statistics 21:619–640. [doi:10.1007/s13253-016-0265-2](https://doi.org/10.1007/s13253-016-0265-2)
- Irvine K.M., Wright W.J., Shanahan E.K., Rodhouse, T.J. (2019) Cohesive framework for modelling plant cover class data. Methods in Ecology and Evolution 10:1749–1760. doi:[10.1111/2041-210X.13262](https://doi.org/10.1111/2041-210X.13262)
- 伊東宏樹 (2020) 植生被度階級データのモデリング. 日本森林学会学術講演集131:S1-5. [doi:10.11519/jfsc.131.0_229](https://doi.org/10.11519/jfsc.131.0_229). コード:[https://github.com/ito4303/jfs131](https://github.com/ito4303/jfs131)
- Itô H. (2020) State-space modeling of the dynamics of temporal plant cover using visually determined class data. PeerJ 8:e9383. [doi:10.7717/peerj.9383](https://doi.org/10.7717/peerj.9383)
- 松浦健太郎 (2016) [StanとRでベイズ統計モデリング](https://www.kyoritsu-pub.co.jp/bookdetail/9784320112421). 共立出版.
- Stan Development Team. (2020) Stan Modeling Language Users Guide and Reference Manual, 2.25.0. [https://mc-stan.org/](https://mc-stan.org/)
