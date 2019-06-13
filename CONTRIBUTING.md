# 开发者须知(CONTRIBUTION GUIDE)
------

首先，非常感谢你能帮助我们完善此项目。但是...，我们都知道，在一个项目的整个
源代码中保持风格和样式的一致性对于项目代码的可读性，可维护性和可用性而言是
非常重要的。所以，我们强烈建议你在为此项目添加代码时能遵循以下风格。

Firstly, I want to say thank you very much for helping me to improve this
project. But..., we all admit that a unified genra and style of the source code
among the whole project is very important to the readability, maintainability
and usability of the project. So we strongly suggest you can adhere to the 
following rules when coding for this project.

## 1. 主要开发语言(Main developing language)

Fortran 90 是本项目的主要开发语言。

Fortran 90 is the main developing language for this project.

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

## 2. Fortran 风格规范(Fotran Style Rules)

### 2.1 文件名(File name)

所有文件都应该以`.F90`或`.f90`为后缀，前者适用于需要预编译的情况，后者适用于
不需要预编译的情况。

All files should be suffixed with `.F90` or `.f90`, the former applies to 
situations where precompilation is required, while the latter applies not.

主程序应该单独放置在一个文件中。除数据模块之外的其他模块应该单独放置在
一个文件中，并以模块名命名文件名。

The main program should be placed in a separate file. Modules, if they are not
data module, should also be placed in a separate file and named after the 
module name.

### 2.2 分号(Semicolons)

不要在行尾加分号, 也不要用分号将两条语句放在同一行。

Do not terminate your lines with semicolons, and do not use semicolons to put
two statements on the same line.

### 2.3 行长度(Line length)

每行不要超过**80个字符**。如果一行写不下，分到多行写，并用'&'符号连接。

Maximum line length is **80 characters**. Break long line into multi ones, and
use '&' line continuation to continue them.

### 2.4 括号(Parentheses)

不要滥用括号。

Use parentheses sparingly.

It is fine, though not required, to use parentheses around tuples. Do not use
them in return statements or conditional statements unless using parentheses for
implied line continuation or to indicate a tuple.

### 2.5 缩进(Indentation)

用**2个空格**来缩进代码。

Indent your code blocks with **2 spaces**. 

绝对不要用tab。也不要把tab和空格混用。对于行连接的情况，你应该要么垂直对齐换行
的元素，或者使用2空格的悬挂式缩进(这时第一行不应该有参数)。

Never use tab or mix tabs and spaces. In cases of implied line continuation,
you should align wrapped elements either vertically, or using a hanging indent 
of 2 spaces, in which case there should be nothing after the open parenthesis or
bracket on the first line.

### 2.6 空行(Blank Lines)

顶级定义之间空两行，如module, subroutine, function等的定义。
变量声明与子程序实体之间空一行。其余结构之间空一行，如功能明显不同的代码结构
之间（如果真出现这种情况，你应该考虑是否需要把两段代码分到不同的子程序中）。

Two blank lines between top-level definitions, such as the module or subroutine
or function definitions. One blank line between variable declaration and 
subprogram entity. One blank line between other structures, such as code
segments with widely different functions (you should first consider if it's
better to separate them into different subprograms if it really happens).

## 2.7 空格(Whitespace)

按照标准的排版规范来使用标点两边的空格。

Follow standard typographic rules for the use of spaces around punctuation.

括号内不要有空格。

No whitespace inside parentheses, brackets or braces.

```fotran
Yes: call foo(arg1(:), arg2(:,:))
```

```fortran
No:  call foo( arg1( : ), arg2( :, : ))
```

不要在逗号前面加空格。但是应该在行尾的行连接符'&'前面加空格。

No whitespace before a comma. Do use whitespace before the line continuation 
character '&' at the end of a line.

在调用有参函数或变量下标时，函数名和参数之间、变量名与下标之间不要有空格。

No whitespace before the open paren/bracket that starts an argument list or 
indexing.

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

### 2.8 Comments and Docstrings

Be sure to use the right style for module, function, method docstrings and
inline comments.


