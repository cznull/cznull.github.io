#include "stdafx.h"

struct node {
	int t;
	int v;
};

//t 0:var 1:con 2:tem 3:fun 4+:sig

int cmp(const char* s1, const char *s2, int n) {
	for (int i = 0; i < n; i++) {
		if (!s1[i] || s1[i] != s2[i]) {
			return 1;
		}
	}
	return s2[n];
}

int gete(std::vector<double> &num,char *s, node &e) {
	char c;
	int i;
	c = *s;
	if ((c <= '9'&&c >= '0') || c == '.') {
		double x = 0, cur = 0.1;
		int type = 0;
		i = 0;
		for (;;) {
			c = s[i];
			if (c <= '9'&&c >= '0') {
				if (type) {
					x += cur * (c - '0');
					cur *= 0.1;
				}
				else {
					x = x * 10.0 + c - '0';
				}
			}
			else if (c == '.') {
				if (type) {
					return -1;
				}
				type = 1;
			}
			else {
				num.push_back(x);
				e.t = 1;
				e.v = num.size() + 1;
				return i;
			}
			i++;
		}
	}
	else if (c <= 'z'&&c >= 'a') {
		i = 0;
		std::vector<char> idf;
		for (;;) {
			c = s[i];
			if (c <= 'z'&&c >= 'a') {
				idf.push_back(c);
			}
			else {
				if (!cmp(idf.data(), "u", idf.size())) {
					e.t = 0;
					e.v = 0;
					return i;
				}
				if (!cmp(idf.data(), "v", idf.size())) {
					e.t = 0;
					e.v = 1;
					return i;
				}
				if (!cmp(idf.data(), "cos", idf.size()) && c == '(') {
					e.t = '(';
					e.v = 1;
					return i + 1;
				}
				if (!cmp(idf.data(), "sin", idf.size()) && c == '(') {
					e.t = '(';
					e.v = 2;
					return i + 1;
				}
				if (!cmp(idf.data(), "sqrt", idf.size()) && c == '(') {
					e.t = '(';
					e.v = 3;
					return i + 1;
				}
				if (!cmp(idf.data(), "ln", idf.size()) && c == '(') {
					e.t = '(';
					e.v = 4;
					return i + 1;
				}
				if (!cmp(idf.data(), "tan", idf.size()) && c == '(') {
					e.t = '(';
					e.v = 5;
					return i + 1;
				}
				if (!cmp(idf.data(), "exp", idf.size()) && c == '(') {
					e.t = '(';
					e.v = 6;
					return i + 1;
				}
				return -1;
			}
			i++;
		}
	}
	else {
		switch (c) {
		case '+':case '-':case '*':case '/':case '(':case ')':
			e.t = c;
			e.v = 0;
			return 1;
		}
	}
}

int tr(std::vector<node> &e, std::vector<double> &num, char *s) {
	int count=0, all=0;
	node temp;
	while (s[all]) {
		count = gete(num, s + all, temp);
		if (count > 0) {
			e.push_back(temp);
			all += count;
		}
		else {
			return -1;
		}
	}
	return 0;
}

int cal(std::vector<node> &e, std::vector<double> &num, oper **ox, int *nox, double **cx, int *outx) {
	node **st;
	int top, i, c, nextresult;
	std::vector<oper> op;
	oper opt;
	c = e.size();
	st = (node**)malloc(c * sizeof(node*));
	top = 0;
	nextresult = num.size() + 2;
	for (i = 0; i < c;) {
		switch (e[i].t) {
		case 0:	case 1: case 2: case '(':
			st[top] = e.data() + i;
			top++;
			i++;
			break;
		case '+':	case '-':
			if (top < 2) {
				st[top] = e.data() + i;
				top++;
				i++;
			}
			else {
				switch (st[top - 1]->t) {
				case 0: case 1: case 2:
					switch (st[top - 2]->t) {
					case '+':
						if (top >= 3) {
							switch (st[top - 3]->t) {
							case 0: case 1: case 2:
								opt.com = st[top - 2]->t;
								opt.num1 = st[top - 3]->v;
								opt.num2 = st[top - 1]->v;
								opt.num3 = nextresult;
								st[top - 3]->t = 2;
								st[top - 3]->v = nextresult;
								nextresult++;
								op.push_back(opt);
								top -= 2;
								break;
							case '(': case '+':	case '-': case '*':	case '/':
								*st[top - 2] = *st[top - 1];
								top--;
								break;
							default:
								goto error;
							}
						}
						else {
							*st[top - 2] = *st[top - 1];
							top--;
						}
						break;
					case '-':
						if (top >= 3) {
							switch (st[top - 3]->t) {
							case 0: case 1: case 2:
								opt.com = st[top - 2]->t;
								opt.num1 = st[top - 3]->v;
								opt.num2 = st[top - 1]->v;
								opt.num3 = nextresult;
								st[top - 3]->t = 2;
								st[top - 3]->v = nextresult;
								nextresult++;
								op.push_back(opt);
								top -= 2;
								break;
							case '(': case '+':	case '-': case '*':	case '/':
								opt.com = 44;
								opt.num1 = st[top - 1]->v;
								opt.num3 = nextresult;
								st[top - 2]->t = 2;
								st[top - 2]->v = nextresult;
								nextresult++;
								op.push_back(opt);
								top--;
								break;
							default:
								goto error;
							}
						}
						else {
							opt.com = 44;
							opt.num1 = st[top - 1]->v;
							opt.num3 = nextresult;
							st[top - 2]->t = 2;
							st[top - 2]->v = nextresult;
							nextresult++;
							op.push_back(opt);
							top--;
							break;
						}
						break;
					case '*':	case '/':
						if (top >= 3) {
							switch (st[top - 3]->t) {
							case 0: case 1: case 2:
								opt.com = st[top - 2]->t;
								opt.num1 = st[top - 3]->v;
								opt.num2 = st[top - 1]->v;
								opt.num3 = nextresult;
								st[top - 3]->t = 2;
								st[top - 3]->v = nextresult;
								nextresult++;
								op.push_back(opt);
								top -= 2;
								break;
							default:
								goto error;
							}
						}
						else {
							goto error;
						}
						break;
					case '(':
						st[top] = e.data() + i;
						top++;
						i++;
						break;
					default:
						goto error;
					}
					break;
				case '(':
					st[top] = e.data() + i;
					top++;
					i++;
					break;
				default:
					goto error;
				}
			}
			break;
		case '*':	case '/':
			if (top < 2) {
				st[top] = e.data() + i;
				top++;
				i++;
			}
			else {
				switch (st[top - 1]->t) {
				case 0: case 1: case 2:
					switch (st[top - 2]->t) {
					case '*':	case '/':
						switch (st[top - 3]->t) {
						case 0: case 1: case 2:
							opt.com = st[top - 2]->t;
							opt.num1 = st[top - 3]->v;
							opt.num2 = st[top - 1]->v;
							opt.num3 = nextresult;
							st[top - 3]->t = 2;
							st[top - 3]->v = nextresult;
							nextresult++;
							op.push_back(opt);
							top -= 2;
							break;
						default:
							goto error;
						}
						break;
					case '(':	case '+':	case '-':
						st[top] = e.data() + i;
						top++;
						i++;
						break;
					}
					break;
				default:
					goto error;
				}
			}
			break;
		case ')':
			while (top >= 2) {
				if (st[top - 2]->t != '(') {
					switch (st[top - 1]->t) {
					case 0: case 1: case 2:
						switch (st[top - 2]->t) {
						case '*': case '/':
							switch (st[top - 3]->t) {
							case 0: case 1: case 2:
								opt.com = st[top - 2]->t;
								opt.num1 = st[top - 3]->v;
								opt.num2 = st[top - 1]->v;
								opt.num3 = nextresult;
								st[top - 3]->t = 2;
								st[top - 3]->v = nextresult;
								nextresult++;
								op.push_back(opt);
								top -= 2;
								break;
							default:
								goto error;
							}
							break;
						case '+':
							if (top >= 3) {
								switch (st[top - 3]->t) {
								case 0: case 1: case 2:
									opt.com = st[top - 2]->t;
									opt.num1 = st[top - 3]->v;
									opt.num2 = st[top - 1]->v;
									opt.num3 = nextresult;
									st[top - 3]->t = 2;
									st[top - 3]->v = nextresult;
									nextresult++;
									op.push_back(opt);
									top -= 2;
									break;
								case '(': case '+':	case '-': case '*':	case '/':
									*st[top - 2] = *st[top - 1];
									top--;
									break;
								default:
									goto error;
								}
							}
							else {
								*st[top - 2] = *st[top - 1];
								top--;
							}
							break;
						case '-':
							if (top >= 3) {
								switch (st[top - 3]->t) {
								case 0: case 1: case 2:
									opt.com = st[top - 2]->t;
									opt.num1 = st[top - 3]->v;
									opt.num2 = st[top - 1]->v;
									opt.num3 = nextresult;
									st[top - 3]->t = 2;
									st[top - 3]->v = nextresult;
									nextresult++;
									op.push_back(opt);
									top -= 2;
									break;
								case '(': case '+':	case '-': case '*':	case '/':
									opt.com = 44;
									opt.num1 = st[top - 1]->v;
									opt.num3 = nextresult;
									st[top - 2]->t = 2;
									st[top - 2]->v = nextresult;
									nextresult++;
									op.push_back(opt);
									top--;
									break;
								default:
									goto error;
								}
							}
							else {
								opt.com = 44;
								opt.num1 = st[top - 1]->v;
								opt.num3 = nextresult;
								st[top - 2]->t = 2;
								st[top - 2]->v = nextresult;
								nextresult++;
								op.push_back(opt);
								top--;
								break;
							}
							break;
						default:
							goto error;
						}
						break;
					default:
						goto error;
					}
				}
				else {
					switch (st[top - 2]->v) {
					case 0:
						top--;
						st[top - 1] = st[top];
						break;
					default:
						opt.com = 256 + st[top - 2]->v;
						opt.num1 = st[top - 1]->v;
						opt.num3 = nextresult;
						st[top - 2]->t = 2;
						st[top - 2]->v = nextresult;
						nextresult++;
						op.push_back(opt);
						top--;
						break;
					}
					i++;
					break;
				}
			}
			break;
		}
	}
	while (top > 1) {
		switch (st[top - 2]->t) {
		case '*': case '/':
			switch (st[top - 3]->t) {
			case 0: case 1: case 2:
				opt.com = st[top - 2]->t;
				opt.num1 = st[top - 3]->v;
				opt.num2 = st[top - 1]->v;
				opt.num3 = nextresult;
				st[top - 3]->t = 2;
				st[top - 3]->v = nextresult;
				nextresult++;
				op.push_back(opt);
				top -= 2;
				break;
			default:
				goto error;
			}
			break;
		case '+':
			if (top >= 3) {
				switch (st[top - 3]->t) {
				case 0: case 1: case 2:
					opt.com = st[top - 2]->t;
					opt.num1 = st[top - 3]->v;
					opt.num2 = st[top - 1]->v;
					opt.num3 = nextresult;
					st[top - 3]->t = 2;
					st[top - 3]->v = nextresult;
					nextresult++;
					op.push_back(opt);
					top -= 2;
					break;
				case '(': case '+':	case '-': case '*':	case '/':
					*st[top - 2] = *st[top - 1];
					top--;
					break;
				}
			}
			else {
				*st[top - 2] = *st[top - 1];
				top--;
			}
			break;
		case '-':
			if (top >= 3) {
				switch (st[top - 3]->t) {
				case 0: case 1: case 2:
					opt.com = st[top - 2]->t;
					opt.num1 = st[top - 3]->v;
					opt.num2 = st[top - 1]->v;
					opt.num3 = nextresult;
					st[top - 3]->t = 2;
					st[top - 3]->v = nextresult;
					nextresult++;
					op.push_back(opt);
					top -= 2;
					break;
				case '(': case '+':	case '-': case '*':	case '/':
					opt.com = 44;
					opt.num1 = st[top - 1]->v;
					opt.num3 = nextresult;
					st[top - 2]->t = 2;
					st[top - 2]->v = nextresult;
					nextresult++;
					op.push_back(opt);
					top--;
					break;
				}
			}
			else {
				opt.com = 44;
				opt.num1 = st[top - 1]->v;
				opt.num3 = nextresult;
				st[top - 2]->t = 2;
				st[top - 2]->v = nextresult;
				nextresult++;
				op.push_back(opt);
				top--;
				break;
			}
			break;
		}
	}
	if (!top) {
		goto error;
	}
	*cx = (double*)malloc(nextresult * sizeof(double));
	memcpy(*cx + 2, num.data(), num.size() * sizeof(double));
	*ox = (oper*)malloc(op.size() * sizeof(oper));
	memcpy(*ox, op.data(), op.size() * sizeof(oper));
	*outx = e[0].v;
	*nox = op.size();
	free(st);
	return 0;
error:
	free(st);
	return -1;
}

int getcoms(express &com, char *x, char *y, char *z) {
	int i;
	std::vector<node> e;
	std::vector<double> num;
	if (tr(e, num, x)) {
		return -1;
	}
	if (cal(e, num, &(com.ox), &(com.nox), &(com.cx), &(com.outx))) {
		return -1;
	}
	e.clear();
	num.clear();
	if (tr(e, num, y)) {
		return -1;
	}
	if (cal(e, num, &(com.oy), &(com.noy), &(com.cy), &(com.outy))) {
		return -1;
	}
	e.clear();
	num.clear();
	if (tr(e, num, z)) {
		return -1;
	}
	if (cal(e, num, &(com.oz), &(com.noz), &(com.cz), &(com.outz))) {
		return -1;
	}
	return 0;
}