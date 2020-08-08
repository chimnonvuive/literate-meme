\LoadClass[12pt, a4paper]{article}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%% Include packages %%%%%%%%%

% fonts and inputs
\usepackage{amsmath, amsfonts, amssymb, gensymb}
\usepackage{mathtools}
\usepackage{times}
%\usepackage{titlesec}
\usepackage[english]{babel}
\usepackage[utf8]{inputenc}
\usepackage{graphicx}
\usepackage{booktabs,multirow, tabulary, float}
\usepackage{array}
\usepackage{longtable}


\usepackage{tikz}
%\usetikzlibrary{arrows.meta} % customizing arrow sizes
%\usetikzlibrary{patterns} % hatching
%\usetikzlibrary{decorations.pathmorphing,patterns} % drawing springs
\usepackage{graphicx}
% graphs in latex: /https://www.bu.edu/math/files/2013/08/tikzpgfmanual.pdf
\usepackage{float}
% for more insights about the float package, check these links: https://tex.stackexchange.com/questions/2275/keeping-tables-figures-close-to-where-they-are-mentioned

% page settings
\usepackage{enumitem}
\usepackage{verbatim}
\usepackage[margin=1in]{geometry}
\usepackage{fancyhdr}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%% Small tweaks %%%%%%%%%%%
\pagestyle{fancy}
%\renewcommand{\sectionmark}[1]{\markright{#1}}
%\renewcommand{\sectionmark}[1]{%
%	\markboth{\thesection\quad #1}{}}
%\fancyhead[RO,RE]{\sectionmark}
%\fancyhf{}
%\chead{\rightmark}
\rhead{\fancyplain{}{\leftmark}}
%\lhead{Engineering Mechanics}
%\rfoot{Page \thepage}

%\renewcommand{\sectionmark}[1]{\markright{#1}}

\DeclarePairedDelimiter\abs{\lvert}{\rvert}
\DeclarePairedDelimiter\norm{\lVert}{\rVert}%

%%%%%% From https://tex.stackexchange.com/questions/43008/absolute-value-symbols %%%%%%
% Swap the definition of \abs* and \norm*, so that \abs and \norm resizes the size of the brackets, and the starred version does not.
%\makeatletter
\let\oldabs\abs
\def\abs{\@ifstar{\oldabs}{\oldabs*}}
\let\oldnorm\norm
\def\norm{\@ifstar{\oldnorm}{\oldnorm*}}
\makeatother

\fancyhead[R]{\includegraphics[width=20pt]{photos/logo.png}}
\renewcommand{\headrulewidth}{0pt}
\renewcommand{\footrulewidth}{0pt}
%\renewcommand{\thesection}{}
%\renewcommand{\thesubsection}{\arabic{section}.\arabic{subsection}}
%\makeatletter
%\def\@seccntformat#1{\csname #1ignore\expandafter\endcsname\csname the#1\endcsname\quad}
%\let\sectionignore\@gobbletwo
%\let\latex@numberline\numberline
%\def\numberline#1{\if\relax#1\relax\else\latex@numberline{#1}\fi}
%\makeatother
%\renewcommand\thesubsection{\@arabic\c@section.\@arabic\c@subsection}
%\makeatother


%\titleformat{\chapter}
%{\normalfont\LARGE\bfseries}{\thechapter}{Chapter}{}
%\titlespacing*{\chapter}{0pt}{3.5ex plus 1ex minus .2ex}{2.3ex plus .2ex}
%\renewcommand{\chaptername}{}
%\titleformat{\chapter}[hang]

\makeatletter         
\def\@maketitle{
	\raggedright
	%\includegraphics[width = 40mm]{logo.jpg}\\[8ex]
	\begin{center}
		{\@title }\\[6ex] 
		{\large  \@author}\\[8ex] 
		\@date\\[8ex]
		\includegraphics[width = 40mm]{photos/logo.png}
\end{center}}
\makeatother


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Convenient notations %%%%%%%




\newcommand{\CC}{\mathbb{C}}
\newcommand{\RR}{\mathbb{R}}
\newcommand{\QQ}{\mathbb{Q}}
\newcommand{\ZZ}{\mathbb{Z}}
\newcommand{\NN}{\mathbb{N}}

\newcommand{\ihat}{\hat{\textbf{i}}}
\newcommand{\jhat}{\hat{\textbf{j}}}
\newcommand{\khat}{\hat{\textbf{k}}}