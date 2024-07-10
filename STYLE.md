# Unofficial PHASTA_NCSU style guide #
(Borrowed from SCOREC CORE style guide. trimmed for needs in PHASTA projects)

The really important rules are:

1. Be clean. If you don't know what clean is, read code until you do.
2. Be consistent. Match the style of surrounding code, unless its not clean.
3. Be concise. Reduce characters and lines unless thats not consistent or clean.

## Files ##

Fortran source has extension `.f`, C++ header files use extension `.h`, C++ source has extension `.cc`,
C source has extension `.c`. 

Files should not exceed 1000 lines,
and lines should not exceed 72 characters.

## Symbols ##

Global symbols are a concatenation of English words.
Ideally, spell the words out fully.
If the word is too long, find a short synonym.
If that doesn't work, find a good abbreviation with some vowels.
For example, `ctor` is better for `constructor` than `cnstrctr`.

Words are concatenated either with nothing between them, in which
case we use camel case to distinguish between words: `createBetterMesh`,
or using underscores, in which case all letters are upper or lower case:
`create_better_mesh`.

All preprocessor symbols and enums should be all upper case,
which forces them to use
underscores to distinguish words: `CREATE_BETTER_MESH`.
Leading or trailing underscores on preprocessor symbols
are reserved for system libraries,
and mortals are forbidden from using these in their symbols:
`__TURN_BACK_HUMAN__`.

Choose one of the above styles, do not mix: `thisIs_Reallynot_good`.

Symbols which are local to a function or structure
don't have to be fully spelled-out words,
since functions and structures should be small,
then there are few symbols to think about,
so things like `n` and `i` are perfectly understandable.

## Braces ##

Either

	for (i = 0; i < 3; ++i) {
	  b[i] = a[i];
	}

or

	for (i = 0; i < 3; ++i)
	{
	  b[i] = a[i];
	}

Anything else, like the GNU style, is too weird.

## Whitespace ##

Currently, most of the code uses 2 spaces for indentation.
Consistency requires that this is followed in existing code,
but new projects may choose a different self-consistent style.

Setup your editor to automate indentation with the style of surrounding code.

Do not mix tabs and spaces, be consistent.

Do not leave trailing whitespace.

## Comments ##

Reserve the use of comments for describing functionality, why some code exists, etc. .

Do not leave commented code in the develop or master branches.  

## Fortran ##

subroutine names should be indicative. The camel case is preferred.
When composing new subroutines, comments of the subroutine functionalities are 
encouraged following the decleration of the subroutine arguments. 

An example

        	subroutine BubASSY()
	!----------------------------------------------------------------------
	!       Called in elmgmr.f
	!       This subroutine is used to assembly bubbles' information after
	!       the loop over blocks on each processor
	!
	!----------------------------------------------------------------------


## C++ ##

Types, including classes, are camel case starting with a capital letter.
Functions and variables should be camel case starting with a lower case letter.
Functions preferrably start with a verb.

A project should be contained within a namespace, that is, all its public
symbols should be in a namespace that is the project prefix.

	namespace apf {

	class SomeType;

	void doSomethingGreat(int a, int b);

	}

## C ##

public symbols should be lowercase with underscores,
starting with the project prefix.

	struct gmi_model;

	double gmi_compute(struct gmi_model* m);

