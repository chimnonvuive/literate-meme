\LoadClass[fleqn,12pt,a4paper,oneside]{report}
\usepackage[utf8]{inputenc}
%\renewcommand{\rmdefault}{unc}
%\usepackage[T1]{fontenc}
\usepackage[english]{babel}
\usepackage{array}
\usepackage{tgtermes}
\usepackage[fleqn]{amsmath}
\usepackage{amsfonts}
\usepackage{amssymb}
\usepackage{biblatex}
\usepackage{booktabs}
%\usepackage{breqn}
\usepackage{caption}
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
\usepackage{nomentbl}
\usepackage{placeins}
\usepackage{pdflscape}
\usepackage[nodisplayskipstretch]{setspace}
\usepackage{siunitx}
\usepackage{stackengine}
\usepackage{subcaption}
\usepackage[table,xcdraw]{xcolor}

%\usepackage{showframe}
\usepackage[tiny]{titlesec}
\renewcommand{\thesection}{\Roman{section}}
\renewcommand{\thesubsection}{\arabic{section}.\arabic{subsection}}
%\titleformat{\chapter}[display]
%{\normalfont%
%	\huge% %change this size to your needs for the first line
%	\bfseries}{\chaptertitlename\ \thechapter}{12pt}{%
%	\Huge %change this size to your needs for the second line
%}


\captionsetup{labelsep=space,justification=raggedright,labelfont=bf,singlelinecheck=false}

\addbibresource{ref.bib}

\newcolumntype{L}[1]{>{\raggedright\let\newline\\\arraybackslash\hspace{0pt}}m{#1}}
\newcolumntype{C}[1]{>{\centering\let\newline\\\arraybackslash\hspace{0pt}}m{#1}}
\newcolumntype{R}[1]{>{\raggedleft\let\newline\\\arraybackslash\hspace{0pt}}m{#1}}

\fancyhf{}
\pagestyle{fancy}


% increase line height
%\renewcommand{\baselinestretch}{1.5}
%\onehalfspacing
\setlength{\belowcaptionskip}{0.6\baselineskip}
\renewcommand{\sectionmark}[1]{\markright{\thesection~ #1}}
\renewcommand{\chaptermark}[1]{\markboth{\thechapter~ #1}{}}


% set default width for figures
\setkeys{Gin}{width=0.6\linewidth}
% set default path for figures
\graphicspath{{./images/}}

\title{Design Project Report}
\author{Group 3}

\makeatletter
\let\thetitle\@title
\let\theauthor\@author
\let\thedate\@date
\makeatother

\newcommand\xrowht[2][0]{\addstackgap[.5\dimexpr#2\relax]{\vphantom{#1}}}
\newcommand{\degc}{^\circ \text{C}}
\newcommand{\s}{\hspace{20mm}}
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

\setlength{\mathindent}{20mm}
\makeglossaries