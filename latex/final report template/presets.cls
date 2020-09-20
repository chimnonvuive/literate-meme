\LoadClass[14pt,a4paper,twoside]{extreport}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{array}
\usepackage[english]{babel}
\usepackage{tgtermes}
\usepackage{amsmath}
\usepackage{amsfonts}
\usepackage{amssymb}
\usepackage[backend=biber]{biblatex}
\usepackage{booktabs}
\usepackage{caption}
\usepackage{fancyhdr}
\usepackage{float}
\usepackage[left=3cm,right=2cm,top=2.5cm,bottom=2.5cm]{geometry}
\usepackage{glossaries}
\usepackage{graphicx}
\usepackage[hidelinks,colorlinks=true,linkcolor=black]{hyperref}
\usepackage{longtable}
\usepackage{multirow}
\usepackage{newtxmath}
\usepackage{newtxtext}
\usepackage{placeins}
\usepackage{setspace}
\usepackage{stackengine}
\usepackage{subcaption}
\usepackage[table,xcdraw]{xcolor}


\renewcommand{\baselinestretch}{1.25}
\newcolumntype{L}[1]{>{\raggedright\let\newline\\\arraybackslash\hspace{0pt}}m{#1}}
\newcolumntype{C}[1]{>{\centering\let\newline\\\arraybackslash\hspace{0pt}}m{#1}}
\newcolumntype{R}[1]{>{\raggedleft\let\newline\\\arraybackslash\hspace{0pt}}m{#1}}

\addbibresource{ref.bib}
\graphicspath{{./images/}}

% set default width for figures
\setkeys{Gin}{width=0.6\linewidth}

%\usepackage{vmargin}
%\setmarginsrb{3 cm}{2.4 cm}{3 cm}{2.4 cm}{1 cm}{1.5 cm}{1 cm}{1.5 cm}

\title{Air Compressor Modeling using MATLAB}
\author{Group 3}

\makeatletter
\let\thetitle\@title
\let\theauthor\@author
\let\thedate\@date
\makeatother

\newcommand\xrowht[2][0]{\addstackgap[.5\dimexpr#2\relax]{\vphantom{#1}}}
\newcommand{\unit}[1]{\ensuremath{\, \mathrm{#1}}}
\newcommand{\unitp}[1]{\ensuremath{\, (\mathrm{#1})}}
\newcommand{\md}[1]{\multicolumn{1}{c}{#1}}
\newcommand{\ih}{\hat{\mathbf{i}}}
\newcommand{\jh}{\hat{\mathbf{j}}}
\newcommand{\kh}{\hat{\mathbf{k}}}
\newcommand{\matgen}[3]{\left[
	\begin{array}{c}
		#1\\
		#2\\
		#3
	\end{array}\right]}

%\pagestyle{fancy}
%\fancyhf{}
%
%\renewcommand{\sectionmark}[1]{\markright{\thesection~ #1}}
%\renewcommand{\chaptermark}[1]{\markboth{\thechapter~ #1}{}}
%
%\fancyhead[LO]{\rightmark}
%\fancyhead[RE]{\leftmark}
%\fancyhead[LE,RO]{\thepage}
%\cfoot{}	
%\renewcommand{\headrulewidth}{0pt}
%%\rhead{\theauthor}
%%\lhead{\thetitle}
%%\cfoot{\thepage}

%\makeglossaries
%\newglossaryentry{latex}
%{
%	name={latex},
%	description={Is a mark up language specially suited 
%		for scientific documents}
%}
%\newglossaryentry{maths}
%{
%	name={mathematics},
%	description={Mathematics is what mathematicians do}
%}
%\newglossaryentry{formula}{
%	name={Chezy equation},
%	description={Chezy equation}
%}
