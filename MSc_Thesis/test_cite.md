---
title: "Untitled"
output:
  html_document:
    keep_md: true
bibliography: grateful-refs.bib
csl: vancouver-superscript.csl
---

@bookdown2016

---
# This is a YAML metadata block, not a front matter
# Cite materials that are cited by flextable like below
# Note that string must be quoted
nocite: '@rmarkdown2018 @rmarkdown2020'
---


``` r
library(ftExtra)
```

```
## Warning: package 'ftExtra' was built under R version 4.4.3
```

```
## Registered S3 method overwritten by 'ftExtra':
##   method                  from     
##   as_flextable.data.frame flextable
```

``` r
data.frame(pkg = 'bookdown @rmarkdown2018 @rmarkdown2020') %>%
  as_flextable() %>%
  colformat_md(pandoc_args = c('--csl', 'harvard-cite-them-right.csl'))
```

```
## Warning in .Deprecated("flextable:::as_flextable.data.frame", msg =
## paste("ftExtra:::as_flextable.data.frame is deprecated", :
## ftExtra:::as_flextable.data.frame is deprecated and will be removed in the
## future release. Consider using flextalbe's implementation by running
## `.S3method("as_flextable", "data.frame", flextable:::as_flextable.data.frame)`
```

---
ftExtra-cite-unnamed-chunk-1: "@rmarkdown2018 @rmarkdown2020"
---

```{=html}
<div class="tabwid"><style>.cl-533491da{}.cl-532b0d36{font-family:'Arial';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-532ef914{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-532fa88c{width:0.75in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 1.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-532fa896{width:0.75in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table data-quarto-disable-processing='true' class='cl-533491da'><thead><tr style="overflow-wrap:break-word;"><th class="cl-532fa88c"><p class="cl-532ef914"><span class="cl-532b0d36">pkg</span></p></th></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-532fa896"><p class="cl-532ef914"><span class="cl-532b0d36">bookdown</span><span class="cl-532b0d36"> </span><span class="cl-532b0d36">Xie,</span><span class="cl-532b0d36"> </span><span class="cl-532b0d36">Allaire</span><span class="cl-532b0d36"> </span><span class="cl-532b0d36">and</span><span class="cl-532b0d36"> </span><span class="cl-532b0d36">Grolemund</span><span class="cl-532b0d36"> </span><span class="cl-532b0d36">(2018)</span><span class="cl-532b0d36"> </span><span class="cl-532b0d36">Xie,</span><span class="cl-532b0d36"> </span><span class="cl-532b0d36">Dervieux</span><span class="cl-532b0d36"> </span><span class="cl-532b0d36">and</span><span class="cl-532b0d36"> </span><span class="cl-532b0d36">Riederer</span><span class="cl-532b0d36"> </span><span class="cl-532b0d36">(2020)</span></p></td></tr></tbody></table></div>
```
