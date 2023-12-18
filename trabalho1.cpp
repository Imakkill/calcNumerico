#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

double funcao(double x) {
    return exp(-pow(x, 2)) - cos(x);
}

double funcaoDerivavada(double x) {
    return sin(x) - 2 * exp(-pow(x, 2));
}

float fi(float x) {
    return pow(x, 3) / 9.0 + 1.0 / 3.0;
}

void executarMetodoBisseccao(ofstream &outputfile, double a, double b, double prec, int n) {
    int k = 0;
    double raiz;
    double x, finicio, meio, fmeio;

    if (fabs(b - a) < prec) {
        x = a;
    } else {
        while (fabs(b - a) > prec && k < n) {
            k++;
            finicio = funcao(a);
            meio = (a + b) / 2;
            fmeio = funcao(meio);
            outputfile << "interação: " << k << ":\t a =" << "\t b= " << b << "\t meio= " << meio << "\t k = " << k << endl;
            if (finicio * fmeio < 0) {
                b = meio;
                raiz = b;
            } else {
                a = meio;
                raiz = a;
            }
        }
    }
    outputfile << "raiz: " << raiz << endl;
    outputfile << "num de int: " << k << endl;
}

void executarMetodoMIL(ofstream &outputfile, double x0, double precisao, int it) {
    int k = 1;
    double raiz, x1, e;

    if (fabs(funcao(x0)) < precisao) {
        e = x0;
        outputfile << "x: " << e << endl;
    } else {
        while (k < it) {
            x1 = fi(x0);
            if (fabs(funcao(x1)) < precisao || fabs(x1 - x0) < precisao || k > it) {
                e = x1;
                outputfile << "x" << k << ": " << e << endl;
                break;
            } else {
                x0 = x1;
                outputfile << "x" << k << ": " << x0 << endl;
                k = k + 1;
            }
        }
    }
    outputfile << "numero de interações: " << k << endl;
}

void executarMetodoNewton(ofstream &outputfile, double x0, double precisao, int it) {
    int k = 1;
    double raiz, fx, fxlinha;

    fx = funcao(x0);
    if (fabs(fx) > precisao) {
        fxlinha = funcaoDerivavada(x0);
        double x1 = x0 - (fx / fxlinha);
        while (fabs(fx) > precisao && fabs(x1 - x0) > precisao && k <= it) {
            k = k + 1;
            x0 = x1;
            fxlinha = funcaoDerivavada(x0);
            x1 = x0 - (fx / fxlinha);
            fx = funcao(x1);
        }
        raiz = x1;
    } else {
        raiz = x0;
    }
    outputfile << "raiz: " << raiz << endl;
    outputfile << "num de it: " << k << endl;
}

void executarMetodoSecante(ofstream &outputfile, double x0, double x1, double prec, int it) {
    int k = 0;
    double raiz, x2;

    if (fabs(funcao(x0)) < prec) {
        raiz = x0;
    } else if (fabs(funcao(x1)) < prec || fabs(x1 - x0) < prec) {
        raiz = x1;
    } else {
        while (fabs(funcao(x1)) > prec && fabs(x1 - x0) > prec && k < it) {
            x2 = x1 - (funcao(x1) * (x1 - x0)) / (funcao(x1) - funcao(x0));
            if (fabs(funcao(x2)) < prec || fabs(x1 - x0) < prec || k >= it) {
                raiz = x2;
                break;
            }
            x0 = x1;
            x1 = x2;
            k = k + 1;
        }
        outputfile << "maior coeficiente: " << raiz << endl;
        outputfile << "numero de interacoes: " << k << endl;
    }
}

void executarMetodoRegulaFalsi(ofstream &outputfile, double a, double b, double prec, int n) {
    int k = 0;
    double raiz, x0, x1, x2;

    if (fabs(b - a) < prec) {
        x2 = a;
    } else {
        x0 = a;
        x1 = b;
        while (fabs(funcao(x1)) > prec && k < n) {
            k++;
            x2 = (x0 * funcao(x1) - x1 * funcao(x0)) / (funcao(x1) - funcao(x0));
            outputfile << "interação: " << k << ":\t a =" << "\t b= " << b << "\t x2= " << x2 << "\t k = " << k << endl;
            if (funcao(x0) * funcao(x2) < 0) {
                x1 = x2;
                raiz = x1;
            } else {
                x0 = x2;
                raiz = x0;
            }
        }
    }
    outputfile << "raiz: " << raiz << endl;
    outputfile << "num de int: " << k << endl;
}

int main() {
    double a = 2, b = 3;
    double x0 = 0.5, x1 = 0.7;
    double prec = 0.1, precisao = 0.0005;
    int n = 100, it = 100;
    int metodo;

    cout << "escolha o metodo (1-bisseccao, 2-MIL, 3-newton, 4-secante, 5-regula falsi, 6-todos)";
    cin >> metodo;

    ofstream outputfile("trab.txt");
    if (outputfile.is_open()) {
        double raiz;
        switch (metodo) {
        case 1: // bissecção
            executarMetodoBisseccao(outputfile, a, b, prec, n);
            break;
        case 2: // MIL
            executarMetodoMIL(outputfile, x0, precisao, it);
            break;
        case 3: // newton
            executarMetodoNewton(outputfile, x0, precisao, it);
            break;
        case 4: // secante
            executarMetodoSecante(outputfile, x0, x1, precisao, it);
            break;
        case 5: // regula falsi
            executarMetodoRegulaFalsi(outputfile, a, b, prec, n);
            break;
        case 6: // todos
            executarMetodoBisseccao(outputfile, a, b, prec, n);
            outputfile << endl;
            executarMetodoMIL(outputfile, x0, precisao, it);
            outputfile << endl;
            executarMetodoNewton(outputfile, x0, precisao, it);
            outputfile << endl;
            executarMetodoSecante(outputfile, x0, x1, precisao, it);
            outputfile << endl;
            executarMetodoRegulaFalsi(outputfile, a, b, prec, n);
            break;
        default:
            cerr << "Opção inválida." << endl;
            break;
        }
        outputfile.close();
    } else {
        cerr << "Erro no arquivo" << endl;
    }

    return 0;
}
