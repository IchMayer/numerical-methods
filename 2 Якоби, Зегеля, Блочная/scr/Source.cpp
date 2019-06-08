#include "includs.h"
#include "Matrix.cpp"

#ifdef _WIN32
#define CLEAR "cls"
void SetColor(int text, int background)
{
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), (WORD)((background << 4) | text));
}
#else
#define CLEAR "clear"
#endif

int k;
bool t;
Matrix<double> a;
int metod;

struct menu
{
public:
	static void circle()
	{
		t = 0;
		metod = 0;
		do {
			system(CLEAR);
			picturemenu();
			printmenu();
			char c = _getch();
			switch (c)
			{
			case 13:
				switch (k)
				{
				case 0:
				{
					system(CLEAR);
					picturestart();
					a.OpenMatrix();
					do
					{
						if (metod == 2)
						{
							int N;
							wcin >> N;
							a.BlockRelax(N);
							break;

						}
#ifdef _DEBUG
						wcout << L"k = " << a.k << endl;
						wcout << L"\t x = ";
						for (size_t i = 0; i < a.result.size(); i++)
							wcout << a.result[i] << "  ";
						wcout << endl;
#endif
						(metod ? a.JacobiMethod() : a.GaussSeidelMethod());
					} while (a.GetLastError() != 2);

					wcout << "\t w = " << a.w << "\t k = " << a.k << endl;
					wcout << setprecision(15);
					for (size_t i = 0; i < a.result.size(); i++)
						wcout << L"\t" << fixed << a.result[i] << L" | " << scientific << i + 1 - a.result[i] << endl;
					wcout << fixed;
					char c1 = _getch();
					while (!c1)
						c1 = _getch();
					break;
				}
				case 1:
					k = 0;
					do
					{
						system(CLEAR);
						picturesettings();
						wcout << endl << L"\te = " << a.e << L"\t w = " << a.w << endl;
						wcout << L"\tmetod :   " <<( metod ? L"Jacobi method" : L"Gauss–Seidel method" )<< endl;
						wcout << L"\tMaxiter = " << a.maxiter << endl << endl;
						printsettings();
						char c = _getch();
						switch (c)
						{
						case 13:
							switch (k)
							{
							case 0:
								cin >> a.e;
								break;
							case 1:
								cin >> a.w;
								break;
							case 2:
							{
								char c1;
								char st[100];
								cin >> c1;
								if(c1 != '0' && c1 != '1')
									cin >> st;
								if (c1 == 'G' || c1 == '0' || c1 == 'g')
									metod = 0;
								else
									if (c1 == '2')
										metod = 2;
									else
										metod = 1;
								break;
							}
							case 3:
								cin >> a.maxiter;
								break;
							case 4:
								k = 10;
								break;
							}
							break;
						case 72:
							k--;
							if (k < 0)
								k = 0;
							break;
						case 80:
							k++;
							if (k > 4)
								k = 4;
							break;
						default:
							break;
						}
					} while (k != 10);
					k = 1;
					break;
				case 2:
					system(CLEAR);
					pictureexit();
					Sleep(5000);
					return;
				default:
					break;
				}
				break;
			case 72:
				k--;
				if (k < 0)
					k = 0;
				break;
			case 80:
				k++;
				if (k > 2)
					k = 2;
				break;
			default:
				break;
			}
		} while (1);
	}
private:
	static void printsettings()
	{
		printstr(L"\t\t   enter e\0", 0);
		printstr(L"\t\t   enter w\0", 1);
		printstr(L"\t\t enter metod\0", 2);
		printstr(L"\t\tenter maxiter\0", 3);
		printstr(L"\t\t back to menu\0", 4);
	}
	static void picturemenu()
	{
		wcout << L"＼　　ヽ　　　　i　　|　　　　　 /　　　/　 \n　　　＼　 \n　　　　　　　　　　　　　　;' ':;,,　　　　 ,;'':;, \n　　　　　　　　　　　　　;'　　 ':;,.,.,,,;'　　';,\n　　ー　　　　　　　　 ,:'　　　　　　　　 　　:::::､ \n　_＿　　　　　　　　,:' ／ 　 　　　　＼ 　　　::::', \n　　　　　二　　　　:'　 ●　　　　　 ●　 　　 ::::::i. \n　　￣　　　　　　　i　 '''　(__人_)　　'''' 　 :::::i \n　　　　-‐　　　　　:　 　　　　　　　　　　　 ::::::i \n　　　　　　　　　　　`:,､ 　　　　　 　 　 :::::::: / \n　　　　／　　　　　　 ,:'　　　　　　　 : :::::::::｀:､ \n　　　　　　　　　　　 ,:'　　　　　　　　 : : :::::::｀:､\n" << endl;
	}
	static void picturesettings()
	{
		wcout << L"　　     　　　_,.. -──- ､, " << endl;
		wcout << L"　　   　　,　'' 　 　　　  `ヽ. " << endl;
		wcout << L"　　　　 ／/¨7__　　/ 　 　 i　 _厂廴 " << endl;
		wcout << L" 　　 /￣( ノ__/　/{　　　　} ｢　（_冫} " << endl;
		wcout << L"　　／￣l＿// 　/-|　 ,!ﾑ ￣|＿｢ ＼＿_ " << endl;
		wcout << L". イ　 ,  /!_∠_ | / /_⊥_,ﾉ ハ　    イ " << endl;
		wcout << L"　/ ／ /  〃ん心 ﾚ'|／　ｆ,心 Y　i ＼_＿＞　 " << endl;
		wcout << L"∠イ 　/   ﾄ弋_ツ　　 　 弋_ﾂ i　 |　 | ＼ " << endl;
		wcout << L"_／ _ノ|　,i　⊂⊃   '      ⊂⊃ ./　 !､＿ン " << endl;
		wcout << L" ￣　　∨|　,小、　` ‐ ' 　　 /|／|　/ " << endl;
		wcout << L"　      Y　|ﾍ＞ 、 ＿ ,.　イﾚ|　 ﾚ' " << endl;
		wcout << L"　　　 r'.|  |;;;入ﾞ亠―亠' );;;;;! 　|､ " << endl;
		wcout << L"　　　 ,ノ:,:|.　!|く  _￣￣￣__У　ﾉ|:,:,ヽ " << endl;
		wcout << L"　　　(:.:.:.:ﾑ人!ﾍ 　` ´ 　　 厂|ノ:.:.:丿" << endl;

	}
	static void pictureexit()
	{
		wcout << L"\t\t\t　　　　　|"<<endl;
		wcout << L"\t\t\t　　　　　|" <<endl ;
		wcout << L"\t\t\t　　　　　|" <<endl ;
		wcout << L"\t\t\t　　　　　|" <<endl ;
		wcout << L"\t\t\t　　　　　|" <<endl ;
		wcout << L"\t\t\t　　　　　|" <<endl ;
		wcout << L"\t\t\t　／￣￣＼|" << endl;
		wcout << L"\t\t\t＜ ´･ 　 |＼" << endl;
		wcout << L"\t\t\t　|　３　| 丶＼" << endl;
		wcout << L"\t\t\t＜ 、･　　|　　＼" << endl;
		wcout << L"\t\t\t　＼＿＿／∪ _ ∪)" << endl;
		wcout << L"\t\t\t　　　　　 Ｕ Ｕ" << endl;
	}
	static void printmenu()
	{
		printstr(L"\t\t\t   start\0", 0);
		printstr(L"\t\t\t settings\0", 1);
		printstr(L"\t\t\t   exit\0", 2);
	}
	static void picturestart()
	{
		wcout << L"膩I嶮薤篝爰曷樔黎㌢´　　｀ⅷ\n";
		wcout << L"艇艀裲f睚鳫巓襴骸　　　　贒憊" << endl;
		wcout << L"殪幢緻I翰儂樔黎夢'”　 　 ,ｨ傾" << endl;
		wcout << L"盥皋袍i耘蚌紕偸′　　　 雫寬I" << endl;
		wcout << L"悗f篝嚠篩i縒縡齢　　 　 Ⅷ辨f" << endl;
		wcout << L"輯駲f迯瓲i軌帶′　　　　　`守I厖孩" << endl;
		wcout << L"幢儂儼巓襴緲′　 　 　 　 　 `守枢i磬廛" << endl;
		wcout << L"嚠篩I縒縡夢'´　　　 　 　 　 　 　 `守峽f" << endl;
		wcout << L"蚌紕襴緲′　　　　　　　　　　　　　　　‘守畝" << endl;
		wcout << L"f瓲軌揄′　　　　　　　　　　　　　,gf毯綴" << endl;
		wcout << L"鳫襴鑿緲　　　　　　　　　　 　 　 奪寔f厦" << endl;
		wcout << L"絨緲′　　　　　　 　 　 　 　　　　 　 ”'罨悳" << endl;
		wcout << L"巓緲′　　　　　　 　 　 　 　 　 　 綴〟 ”'罨椁" << endl;
		wcout << L"巓登嶮 薤篝㎜㎜ g　 　 緲　 　 甯體i爺綴｡, ”'罨琥" << endl;
		wcout << L"I軌襴暹 甯幗緲fi'　　 緲',纜　　贒i綟碕碚爺綴｡ ”'罨皴" << endl;
		wcout << L"巓襴驫 霤I緲緲　　 纜穐　　甯絛跨飩i髢綴馳爺綴｡`'等誄" << endl;
	}
	static void printstr(const wchar_t a[], int number)
	{
		if (k == number)
		{
			int i = 0;
#ifdef _WIN32
			SetColor(5, 0);
#endif
			while (a[i] == L'\t')
				wcout << L'\t', i++;
			while (a[i] == L' ')
				wcout << L' ', i++;
			for (; i < wcslen(a); i++)
			{
				char c = (int)a[i] - (int)L'a' + (int)L'A';
				wcout << c;
			}
			wcout << endl;
#ifdef _WIN32
			SetColor(15, 0);
#endif
		}
		else
		{
			wcout << a << endl;
		}
	}
};

int wmain(int argc, wchar_t* argv[])
{
#ifdef _WIN32
	SetConsoleTitle("Numecal methods Mayer Valera lab #2");
#else
	cout << "\033]0;" << "Numecal methods Mayer Valera lab #2" << "\007";
#endif 
	_setmode(_fileno(stdout), _O_U16TEXT);
	menu::circle();
	Matrix<double> a;
	return 0;
}