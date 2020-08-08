\LoadClass[14pt,a4paper,oneside]{extreport}
\usepackage{extsizes}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{pdflscape}
\usepackage[russian,english]{babel}
\usepackage{amsmath}
\usepackage{amsfonts}
\usepackage{amssymb}
\usepackage{graphicx}
\usepackage{setspace}
\usepackage{newtxtext,newtxmath,multirow}
\usepackage{stackengine,booktabs,float,placeins,longtable}
\usepackage[margin=2.4cm]{geometry}
\usepackage[table,xcdraw]{xcolor}
\usepackage[none]{hyphenat}
\renewcommand{\baselinestretch}{1.5}
\newcommand\xrowht[2][0]{\addstackgap[.5\dimexpr#2\relax]{\vphantom{#1}}}
\newcommand{\unit}[1]{\ensuremath{\, \mathrm{#1}}}
\newcommand{\md}[1]{\multicolumn{1}{c}{#1}}
\usepackage{array}
\newcolumntype{L}[1]{>{\raggedright\let\newline\\\arraybackslash\hspace{0pt}}m{#1}}
\newcolumntype{C}[1]{>{\centering\let\newline\\\arraybackslash\hspace{0pt}}m{#1}}
\newcolumntype{R}[1]{>{\raggedleft\let\newline\\\arraybackslash\hspace{0pt}}m{#1}}

\graphicspath{{images/}}
\usepackage{parskip}
\usepackage{fancyhdr}
\usepackage{vmargin}
\setmarginsrb{25mm}{20mm}{15mm}{20mm}{20mm}{1.5 cm}{1 cm}{1.5 cm}

\title{Lab Report}
\author{Nguyen Quy Khoi}
\date{21 July 2020}

\makeatletter
\let\thetitle\@title
\let\theauthor\@author
\let\thedate\@date
\makeatother

\pagestyle{fancy}
\fancyhf{}
\rhead{\theauthor}
\lhead{\thetitle}
\cfoot{\thepage}