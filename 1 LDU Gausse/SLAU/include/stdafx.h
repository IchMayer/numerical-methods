// stdafx.h: включаемый файл для стандартных системных включаемых файлов
// или включаемых файлов для конкретного проекта, которые часто используются, но
// не часто изменяются
//

#pragma once

#include "targetver.h"

#define WIN32_LEAN_AND_MEAN             // Исключите редко используемые компоненты из заголовков Windows
// Файлы заголовков Windows:
#include <windows.h>

// Файлы заголовков C RunTime
#include <stdlib.h>
#include <malloc.h>
#include <memory.h>
#include <tchar.h>

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <time.h>
#include <commdlg.h>
#include <comdef.h>

#define max(a,b) (a > b ? a : b)
#define min(a,b) (a > b ? b : a)
#define abs(a) (a < 0 ? -a: a)

using namespace std;
