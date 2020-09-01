\LoadClass{report}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[russian,english]{babel}
\usepackage{amsmath}
\usepackage{amsfonts}
\usepackage{amssymb}
\usepackage{graphicx}
\usepackage{setspace}
\usepackage{newtxtext,newtxmath,multirow}
\usepackage{stackengine,booktabs,float,placeins,longtable}
\usepackage[a5paper,top=25mm,left=20mm,right=20mm,bottom=25mm]{geometry}
\usepackage[table,xcdraw]{xcolor}
\usepackage[none]{hyphenat}

\newcommand\xrowht[2][0]{\addstackgap[.5\dimexpr#2\relax]{\vphantom{#1}}}
\newcommand{\unit}[1]{\ensuremath{\, \mathrm{#1}}}
\newcommand{\unitp}[1]{\ensuremath{\, (\mathrm{#1})}}
\newcommand{\md}[1]{\multicolumn{1}{c}{#1}}
\usepackage{array}
\newcolumntype{L}[1]{>{\raggedright\let\newline\\\arraybackslash\hspace{0pt}}m{#1}}
\newcolumntype{C}[1]{>{\centering\let\newline\\\arraybackslash\hspace{0pt}}m{#1}}
\newcolumntype{R}[1]{>{\raggedleft\let\newline\\\arraybackslash\hspace{0pt}}m{#1}}
\usepackage[hidelinks]{hyperref}
\graphicspath{{images/}}
\usepackage{parskip}
\usepackage{fancyhdr}
\fancyhf{}
\pagestyle{fancy}
\renewcommand{\sectionmark}[1]{\markright{\thesection~ #1}}
\renewcommand{\chaptermark}[1]{\markboth{\thechapter~ #1}{}}

\fancyhead[LO]{\rightmark}
\fancyhead[RE]{\leftmark}
\fancyhead[LE,RO]{\thepage}
\cfoot{}	
\renewcommand{\headrulewidth}{0pt}

\title{Assignment}
\author{Nguyen Quy Khoi}
\date{\today}

\newcommand{\vb}[2]{\hat{\bm{#1}}_{\bm{#2}}}
\newcommand{\ih}{\hat{\mathbf{i}}}
\newcommand{\jh}{\hat{\mathbf{j}}}
\newcommand{\kh}{\hat{\mathbf{k}}}

\makeatletter
\let\thetitle\@title
\let\theauthor\@author
\let\thedate\@date
\makeatother