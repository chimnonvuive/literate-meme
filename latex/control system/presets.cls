\LoadClass[12pt,a4paper,twoside]{report}
\usepackage[utf8]{inputenc}
%\renewcommand{\rmdefault}{unc}
%\usepackage[T1]{fontenc}
\usepackage[english]{babel}
\usepackage{minted}
\usepackage{array}
\usepackage{tgtermes}
\usepackage{amsmath}
\usepackage{amsfonts}
\usepackage{amssymb}
\usepackage{biblatex}
\usepackage{booktabs}
\usepackage{caption}
\usepackage{mhchem}
\usepackage{etoolbox}
\usepackage{fancyhdr}
\usepackage{float}
\usepackage[left=3cm,right=2cm,top=2.5cm,bottom=2.5cm]{geometry}
\usepackage[toc,automake]{glossaries}
\usepackage{graphicx}
\usepackage[hidelinks,colorlinks=true,linkcolor=black,citecolor=black]{hyperref}
\usepackage{imakeidx}
\usepackage{lipsum}
\usepackage{longtable}
\usepackage{multicol}
\usepackage{multirow}
\usepackage{newtxmath}
\usepackage{newtxtext}
\usepackage{placeins}
\usepackage{setspace}
\usepackage{stackengine}
\usepackage{subcaption}
\usepackage[table,xcdraw]{xcolor}

\addbibresource{ref.bib}
\newcolumntype{L}[1]{>{\raggedright\let\newline\\\arraybackslash\hspace{0pt}}m{#1}}
\newcolumntype{C}[1]{>{\centering\let\newline\\\arraybackslash\hspace{0pt}}m{#1}}
\newcolumntype{R}[1]{>{\raggedleft\let\newline\\\arraybackslash\hspace{0pt}}m{#1}}

\fancyhf{}
\pagestyle{fancy}


% increase line height
%\renewcommand{\baselinestretch}{1.25}
%\setlength{\belowcaptionskip}{0.6\baselineskip}
\renewcommand{\sectionmark}[1]{\markright{\thesection~ #1}}
\renewcommand{\chaptermark}[1]{\markboth{\thechapter~ #1}{}}


% set default width for figures
\setkeys{Gin}{width=0.6\linewidth}
% set default path for figures
\graphicspath{{./images/}}

\title{PID Controller Design for Self-excited Brushed DC Motor}
\author{Nguyen Quy Khoi}

\makeatletter
\let\thetitle\@title
\let\theauthor\@author
\let\thedate\@date
\makeatother

\newcommand\xrowht[2][0]{\addstackgap[.5\dimexpr#2\relax]{\vphantom{#1}}}
\newcommand{\degc}{^\circ C}
\newcommand{\unit}[1]{\ensuremath{\, \mathrm{#1}}}
\newcommand{\unitp}[1]{\ensuremath{\, (\mathrm{#1})}}
\newcommand{\md}[1]{\multicolumn{1}{c}{#1}}
\newcommand{\ih}{\hat{\mathbf{i}}}
\newcommand{\jh}{\hat{\mathbf{j}}}
\newcommand{\kh}{\hat{\mathbf{k}}}
\newcommand{\matgen}[3]{
	\left[
	\begin{array}{c}
		#1\\
		#2\\
		#3
	\end{array}
	\right]}

\fancyhead[LO]{\rightmark}
\fancyhead[RE]{\leftmark}
\fancyhead[LE,RO]{\thepage}
%\cfoot{}	
\renewcommand{\headrulewidth}{0pt}
%\rhead{\theauthor}
%\lhead{\thetitle}
%\cfoot{\thepage}

%\setlength{\mathindent}{20mm}

\makeindex[intoc]
\makeglossaries