\LoadClass[13pt,a4paper,twoside]{extreport}
\usepackage[utf8]{inputenc}
%\usepackage[T1]{fontenc}
%\usepackage[utf8]{vietnam}
\usepackage{array}
\usepackage[english]{babel}
\usepackage{tgtermes}
\usepackage{amsmath}
\usepackage{amsfonts}
\usepackage{amssymb}
\usepackage{booktabs}
\usepackage{float}
\usepackage[left=3cm,right=2cm,top=2.5cm,bottom=2.5cm]{geometry}
\usepackage{graphicx}
\usepackage[hidelinks,colorlinks=true,linkcolor=blue,linktocpage=true]{hyperref}
\usepackage{longtable}
\usepackage{multirow}
\usepackage{newtxmath}
\usepackage{newtxtext}
\usepackage{placeins}
\usepackage{setspace}
\usepackage{stackengine}
\usepackage[nottoc]{tocbibind}
\usepackage[table,xcdraw]{xcolor}

\renewcommand{\baselinestretch}{1.25}
\newcommand\xrowht[2][0]{\addstackgap[.5\dimexpr#2\relax]{\vphantom{#1}}}
\newcommand{\unit}[1]{\ensuremath{\, \mathrm{#1}}}
\newcommand{\unitp}[1]{\ensuremath{\, (\mathrm{#1})}}
\newcommand{\md}[1]{\multicolumn{1}{c}{#1}}
\newcolumntype{L}[1]{>{\raggedright\let\newline\\\arraybackslash\hspace{0pt}}m{#1}}
\newcolumntype{C}[1]{>{\centering\let\newline\\\arraybackslash\hspace{0pt}}m{#1}}
\newcolumntype{R}[1]{>{\raggedleft\let\newline\\\arraybackslash\hspace{0pt}}m{#1}}

%\addbibresource{ref.bib}
\graphicspath{{images/}}
%\usepackage{parskip}
\usepackage{fancyhdr}
%\usepackage{vmargin}
%\setmarginsrb{3 cm}{2.4 cm}{3 cm}{2.4 cm}{1 cm}{1.5 cm}{1 cm}{1.5 cm}

\title{Air Compressor Modeling using MATLAB}
\author{Group 3}
\date{30 June 2020}

\makeatletter
\let\thetitle\@title
\let\theauthor\@author
\let\thedate\@date
\makeatother

\newcommand{\ih}{\hat{\mathbf{i}}}
\newcommand{\jh}{\hat{\mathbf{j}}}
\newcommand{\kh}{\hat{\mathbf{k}}}
\newcommand{\matgen}[3]{\left[
	\begin{array}{c}
		#1\\
		#2\\
		#3
	\end{array}\right]}

\geometry{a4paper}
\pagestyle{fancy}
\fancyhf{}

\renewcommand{\sectionmark}[1]{\markright{\thesection~ #1}}
\renewcommand{\chaptermark}[1]{\markboth{\thechapter~ #1}{}}

\fancyhead[LO]{\rightmark}
\fancyhead[RE]{\leftmark}
\fancyhead[LE,RO]{\thepage}
\cfoot{}	
\renewcommand{\headrulewidth}{0pt}
%\rhead{\theauthor}
%\lhead{\thetitle}
%\cfoot{\thepage}