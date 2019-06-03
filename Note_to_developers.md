# NOTE TO DEVELOPERS
------

Firstly, I want to say thank you very much for helping me to improve this
project. But..., we all admit that a unified genra and style of the source code
among the whole project is very important to the readability, maintainability
and usability of the project. So we strongly suggest you can adhere to the 
following rules when coding.

## 1. Main developing language

Fortran 90 is our main developing language.
Of course Fortran 90 cannot do everything, and it can not solve all problems 
efficiently. So, other languages can be chosen when there is ready-made source
code, such as the specialized math functions. These verified source code in 
other languages can be used in this project without any changes. But if you do 
this, make sure the licenses of them are all compatible with that in this
project and make acknowledgement in your source code file and let us know this
when you make a commitment.

The priority sequence of the source code languages of this project are:
1. Fortran 90
2. C/C++
3. Python

The following rules for the genre and style of this project's source code are
defined for Fortran 90 languages, but you can make an analogy for other
languages, such as Fortran 77 and C/C++. 

Once more, we strongly suggest you can adhere to the following rules when
coding.

## 2. Fotran Style Rules

### 2.1 Semicolons

Do not terminate your lines with semicolons, and do not use semicolons to put
two statements on the same line.

### 2.2 Line length

Maximum line length is *80 characters*. Use '&' line continuation for long lines.

### 2.3 Parentheses

Use parentheses sparingly.

It is fine, though not required, to use parentheses around tuples. Do not use
them in return statements or conditional statements unless using parentheses for
implied line continuation or to indicate a tuple.

### 2.4 Indentation

Indent your code blocks with *2 spaces*. 

Never use tab or mix tabs and spaces. In cases of implied line continuation,
you should align wrapped elements either vertically, or using a hanging indent 
of 2 spaces, in which case there should be nothing after the open parenthesis or
bracket on the first line.

### 2.5 Blank Lines

Two blank lines between top-level definitions, such as the module or subroutine
or function definitions. One blank line between variable declaration and 
line and the first method. No blank line following a `def` line. Use single
blank lines as you judge appropriate within functions or methods.

## 2.6 Whitespace

Follow standard typographic rules for the use of spaces around punctuation.

No whitespace inside parentheses, brackets or braces.

```fotran
Yes: call foo(arg1(:), arg2(:,:))
```

```fortran
No:  call foo( arg1( : ), arg2( :, : ))
```

No whitespace before a comma, semicolon, or colon. Do use whitespace before the
line continuation character '&' at the end of a line.

No whitespace before the open paren/bracket that starts an argument list,
indexing or slicing.

```fortran
Yes: call foo(arg1(1))
```

```fortran
No:  call foo (arg (1))
```

Surround binary operators with a single space on either side for comparisons 
(`==, <, >, /=, <>, <=, >=`), and Booleans (`.and., .or., .not.`). 
Use your better judgment for the insertion of spaces around assignments (`=`)
and arithmetic operators (`+`, `-`, `*`, `/`, `**`).

```fortran
Yes: x == 1
```

```fortran
No:  x<1
```

Don't use spaces to vertically align tokens on adjective lines and on
consecutive lines, since it becomes a maintenance burden (applies to `=` etc.):

```fortran
Yes:
  foo = 1000  ! comment
  long_name = 2  ! comment that should not be aligned
```

```fortran
No:
  foo       = 1000  ! comment
  long_name = 2     ! comment that should not be aligned
```

### 2.7 Comments and Docstrings

Be sure to use the right style for module, function, method docstrings and
inline comments.


