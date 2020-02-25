\LoadClass{beamer}
%
% Choose how your presentation looks.
%
% For more themes, color themes and font themes, see:
% http://deic.uab.es/~iblanes/beamer_gallery/index_by_theme.html
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% better vector notation by Todd Lehman
% visit https://tex.stackexchange.com/questions/44017/dot-notation-for-derivative-of-a-vector for more information
\usepackage[T1]{fontenc}  % For correct {}s in \texttt
\usepackage{amsmath}      % For better \dot placement
\usepackage[b]{esvect}    % For \vv
\usepackage{tikz}         % For arrow and dots in \xvec

% --- Macro \xvec
\makeatletter
\newlength\xvec@height%
\newlength\xvec@depth%
\newlength\xvec@width%
\newcommand{\xvec}[2][]{%
\ifmmode%
\settoheight{\xvec@height}{$#2$}%
\settodepth{\xvec@depth}{$#2$}%
\settowidth{\xvec@width}{$#2$}%
\else%
\settoheight{\xvec@height}{#2}%
\settodepth{\xvec@depth}{#2}%
\settowidth{\xvec@width}{#2}%
\fi%
\def\xvec@arg{#1}%
\def\xvec@dd{:}%
\def\xvec@d{.}%
\raisebox{.2ex}{\raisebox{\xvec@height}{\rlap{%
			\kern.05em%  (Because left edge of drawing is at .05em)
			\begin{tikzpicture}[scale=1]
			\pgfsetroundcap
			\draw (.05em,0)--(\xvec@width-.05em,0);
			\draw (\xvec@width-.05em,0)--(\xvec@width-.15em, .075em);
			\draw (\xvec@width-.05em,0)--(\xvec@width-.15em,-.075em);
			\ifx\xvec@arg\xvec@d%
			\fill(\xvec@width*.45,.5ex) circle (.5pt);%
			\else\ifx\xvec@arg\xvec@dd%
			\fill(\xvec@width*.30,.5ex) circle (.5pt);%
			\fill(\xvec@width*.65,.5ex) circle (.5pt);%
			\fi\fi%
			\end{tikzpicture}%
}}}%
#2%
}
\makeatother

% --- Override \vec with an invocation of \xvec.
\let\stdvec\vec
\renewcommand{\vec}[1]{\xvec[]{#1}}
% --- Define \dvec and \ddvec for dotted and double-dotted vectors.
\newcommand{\dvec}[1]{\xvec[.]{#1}}
\newcommand{\ddvec}[1]{\xvec[:]{#1}}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newcommand{\vb}[2]{\vec{\bm{#1}}_{\bm{#2}}}
\newcommand{\ih}{\vec{\mathbf{i}}}
\newcommand{\jh}{\vec{\mathbf{j}}}
\newcommand{\kh}{\vec{\mathbf{k}}}
\mode<presentation>
{
\usetheme{CambridgeUS}
\usecolortheme{default}
\usefonttheme{professionalfonts}
\setbeamertemplate{navigation symbols}{}
\setbeamertemplate{caption}[numbered]
} 

%\usepackage{helvet}
\usepackage[english]{babel}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{wasysym}
\usepackage{bm, listings, xcolor, matlab-prettifier}

%\renewcommand{\familydefault}{\sfdefault}

\title[Instructor: Dr. Nguyen Tan Tien]{Velocity - Acceleration Analysis}
\author[Khoi]{Nguyen Quy Khoi}
\institute[HCMUT-OISP-FME]{HCMUT-OISP - Faculty of Mechanical Engineering}
\date{20 February, 2020}