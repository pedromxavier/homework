\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{homework}[2020/06/30 Homework LaTeX class]
\LoadClass{article}

\RequirePackage{amsthm}
\RequirePackage{amsmath}
\RequirePackage{ifthen}
\RequirePackage{amssymb}
\RequirePackage{mathrsfs}
\RequirePackage{xcolor}
\RequirePackage{listings}
\RequirePackage[utf8]{inputenc}
\RequirePackage[brazil]{babel}
\RequirePackage{multicol}
\RequirePackage{float}
\RequirePackage{enumitem}
\RequirePackage{textgreek}

%% Pgf - Tikz
\RequirePackage{tikz}
\RequirePackage{pgfplots}

\newcommand{\vertex}{\node[vertex]}

\newcommand{\cqd}{%%
	\ensuremath{\blacksquare}
}

\newcommand{\gvertex}[4][]{%%
	\node (#2) at (#3, #4) {};
	\draw [#1] (#2) circle [radius=2pt];
}

\newcommand{\glink}[3][gray!30]{%%
	\draw [#1] (#2) -- (#3);
}

\setlength{\parindent}{0em}
\setlength{\parskip}{1em}
\setlength{\voffset}{-1in}
\setlength{\textheight}{640pt}

\allowdisplaybreaks

\definecolor{codecolor}{rgb}{0.0, 0.0, 0.6}
\definecolor{stringcolor}{rgb}{0.0, 0.6, 0.0}
\definecolor{emphcolor}{rgb}{0.6, 0.0, 0.6}
\definecolor{commentcolor}{rgb}{0.6, 0.6, 0.6}

\newcommand{\keywordstyle}{%%
	\itshape\bfseries\color{codecolor}
}

\newcommand{\emphstyle}{%%
	\itshape\bfseries\color{emphcolor}
}

\newcommand{\stringstyle}{%%
	\itshape\bfseries\color{stringcolor}
}

\newcommand{\commentstyle}{%%
	\itshape\bfseries\color{commentcolor}
}

\newcommand{\code}[1]{%%
	{\ttfamily\color{codecolor}\bfseries #1}%%
}

\newcommand{\type}[1]{%%
	{\emphstyle #1}%%
}

\renewcommand{\vec}[1]{%%
	\ensuremath{\mathbf{#1}}
}

\renewcommand{\qedsymbol}{%%
	\ensuremath{\blacksquare}
}

\DeclareRobustCommand{\rchi}{{\mathpalette\irchi\relax}}
\newcommand{\irchi}[2]{\raisebox{\depth}{$#1\chi$}} % inner command, used by \rchi

\DeclareRobustCommand{\rsqrt}{{\mathpalette\irsqrt\relax}}
\newcommand{\irsqrt}[2]{\raisebox{\depth}{$#1\sqrt{}$}} % inner command, used by \rchi

\newcounter{algors}
\setcounter{algors}{0}

\lstdefinestyle{algorstyle}
{
	numberblanklines=true,
	tabsize=4,
	gobble=4,
	frame=single,
	stringstyle={\stringstyle},
	keywordstyle={\keywordstyle},
	morekeywords={%  
		def, lambda, seja, sejam, enquanto, para, de, cada, em, se, retorna, continua, e, ou%
	},
	emphstyle={\emphstyle},
	emph={%
	},
	commentstyle={\commentstyle},
	numbers=left,
	stepnumber=1,
	showstringspaces=false,
	basicstyle=\ttfamily,
	mathescape=true,
	breaklines=true,
	breakatwhitespace=false,
	extendedchars=true,
	inputencoding=utf8,
	morecomment=[l]{//},
	morecomment=[s]{/*}{*/},
	morestring=[b]",
	literate=
	{ã}{{\~a}}1                         {õ}{{\~o}}1 
	{Ã}{{\~A}}1                         {Õ}{{\~O}}1 
	{á}{{\'a}}1 {é}{{\'e}}1 {í}{{\'i}}1 {ó}{{\'o}}1 {ú}{{\'u}}1
	{Á}{{\'A}}1 {É}{{\'E}}1 {Í}{{\'I}}1 {Ó}{{\'O}}1 {Ú}{{\'U}}1
	{à}{{\`a}}1 {è}{{\`e}}1 {ì}{{\`i}}1 {ò}{{\`o}}1 {ù}{{\`u}}1
	{À}{{\`A}}1 {È}{{\'E}}1 {Ì}{{\`I}}1 {Ò}{{\`O}}1 {Ù}{{\`U}}1
	{ä}{{\"a}}1 {ë}{{\"e}}1 {ï}{{\"i}}1 {ö}{{\"o}}1 {ü}{{\"u}}1
	{Ä}{{\"A}}1 {Ë}{{\"E}}1 {Ï}{{\"I}}1 {Ö}{{\"O}}1 {Ü}{{\"U}}1
	{â}{{\^a}}1 {ê}{{\^e}}1 {î}{{\^i}}1 {ô}{{\^o}}1 {û}{{\^u}}1
	{Â}{{\^A}}1 {Ê}{{\^E}}1 {Î}{{\^I}}1 {Ô}{{\^O}}1 {Û}{{\^U}}1
	{œ}{{\oe}}1 {Œ}{{\OE}}1 {æ}{{\ae}}1 {Æ}{{\AE}}1 {ß}{{\ss}}1
	{ç}{{\c c}}1 {Ç}{{\c C}}1 {ø}{{\o}}1 {å}{{\r a}}1 {Å}{{\r A}}1
	{€}{{\EUR}}1 {£}{{\pounds}}1
	%%
	{~}{{{\raisebox{0.5ex}{\texttildelow}}}}1
	{^}{{{\raisebox{-0.5ex}{\textasciicircum}}}}1
	{<-}{{{$\leftarrow$}}}1
	{->}{{{$\rightarrow$}}}1
	{<->}{{{$\longleftrightarrow$}}}1
	%%
	{até}{{{\keywordstyle at{\'e}}}}3
	{senão}{{{\keywordstyle sen{\~a}o}}}5
	%%
	{\_fila\_}{{{\emphstyle fila}}}4
	{\_pilha\_}{{{\emphstyle pilha}}}5
	{\_vazio\_}{{{\emphstyle vazio}}}5
	{\_v\_}{{{\emphstyle v}}}1
	{\_e\_}{{{\emphstyle e}}}1
	{\_nulo\_}{{{\emphstyle nulo}}}4
	{\_grau\_}{{{\emphstyle grau}}}4
	{\_soma\_}{{{\emphstyle soma}}}4
	{\_viz\_}{{{\emphstyle viz}}}3
	{\_max\_}{{{\emphstyle max}}}3
	{\_min\_}{{{\emphstyle min}}}3
	{\_mergesort\_}{{{\emphstyle merge-sort}}}{9}
	,
}

\lstnewenvironment{algor}[1]{%%
	\stepcounter{algors}
	\lstset{style=algorstyle, title={Algoritmo \thealgors: #1}}
}{%%
}

\newcommand{\bft}[1]{%%
	\ensuremath{\text{\bfseries#1}}
}

\DeclareMathOperator*{\argmax}{arg\,max}
\DeclareMathOperator*{\argmin}{arg\,min}

%% EOI
\endinput