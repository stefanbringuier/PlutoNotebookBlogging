<!--
Add here global page variables to use throughout your website.
-->
+++
author = "Stefan Bringuier"
mintoclevel = 2


# Add here files or directories that should be ignored by Franklin, otherwise
# these files might be copied and, if markdown, processed by Franklin which
# you might not want. Indicate directories by ending the name with a `/`.
# Base files such as LICENSE.md and README.md are ignored by default.
ignore = ["node_modules/"]

# RSS (the website_{title, descr, url} must be defined to get RSS)
generate_rss =  false
website_title = "Pluto.jl Notebook Blogs - Stefan Bringuier"
website_descr = "Blogging various topics using Pluto.jl"
website_url   = "https://stefanbringuier.github.io/PlutoNotebookBlogging"
+++

@def prepath = "PlutoNotebookBlogging"


<!--
Add here global latex commands to use throughout your pages.
-->
\newcommand{\R}{\mathbb R}
\newcommand{\scal}[1]{\langle #1 \rangle}
\newcommand{\note}[1]{@@note @@title warning Note@@ @@content #1 @@ @@} \newcommand{\warn}[1]{@@warning @@title warning Warning!@@ @@content #1 @@ @@}