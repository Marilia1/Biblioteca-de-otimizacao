#include <iostream>
#include <list>

#ifndef __MDL_H__
#define __MDL_H__

//-----------------Programação Dinamica------------------
//Os algoritmos trabalhados na programação dinamica sera:
//Triângulo de Pascal
//Troco em moedas
//Mochila binária

/*Aluna: Marília Karla Soares Lopes - 11413590
O primeiro programa é o triângulo de Pascal onde criamos uma função 
que recebem como parametro as combinações que desejamos procurar no
triangulo, com ela sabaremos qual numero que pertence aquela linha e coluna */

int min(int a, int b);
/*Função que recebe como parametro a linha e coluna presente no triangulo*/
int Coeficiente(int lin, int col);

/*Aluna: Marília Karla Soares Lopes - 11413590
O segundo programa é o troco em moedas, esse programa nos diz quantas 
possibilidades de combinação para um certo valos de troco podemos ter 
de acordo com os valores de moedas disponíveis */

/*Função que recebe as moedas que temos disponiveis para dar o troco desejado */
int count( int S[], int m, int n );

/*Aluna: Marília Karla Soares Lopes - 11413590
O terceiro programa é o de mochila binária, onde a ideia é que temos uma mchila que 
tem a capacidade para suportar um peso n, e temos objetos com diferentes pesos, e queremos
saber quais itens a mochila aguenta*/

int max(int a, int b) { return (a > b)? a : b; }

// Retrona o valor máximo que pode ser colocado na mochila com uma capacidade C
int knapSack(int C, int pitem[], int val[], int n)



//--------------------Algoritmos gulosos----------------
//Os algoritmos trabalhados nos gulosos serão:
//Mochila fracionária
//Coloração de grafos

/*Aluna: Marília Karla Soares Lopes - 11413590
O primeiro programa é a mochila fracionária, é basicamente a mesma ideia do algoritmo 
mostrado acima, porem dessa vez caso um item tenha um peso que a mochila não suporta
apenas parte dele será add obtendo o valor total da capacidade C da mochila*/

struct Item //estrutura que representa item
{
    int valor, peso;
 
    // Construtor
    Item(int valor, int peso) : valor(valor), peso(peso)
    {}
};

// Função que compara os intes de acordo com a relação do seu valor pelo peso
bool cmp(struct Item a, struct Item b);

return r1 > r2; //retorna resultado da relação entre os pesos e compara r1 com r2
/* Função principal que recebe como parametro a capacidade da mochila e os pesos dos itens
que se deseja colocar dentro da mochila*/
double fractionalKnapsack(int C, struct Item arr[], int n)

// retorna o valor final
return valorfinal;



/*Aluna: Marília Karla Soares Lopes - 11413590
O segundo programa é Coloração de grafos,onde o problema é, dado m cores, 
temos que encontrar uma maneira de colorir os vértices de um grafo de modo que 
não tenha dois vértices adjacentes com a mesma cor. Esse algoritmo não garante 
utilizar cores mínimas, mas garante um limite superior do número de cores. 
O algoritmo básico nunca usa mais do que d cores + 1, onde d é o grau máximo 
de um vértice no gráfico dado. */

//Classe que representa um grafo não direcionavel
class Grafo
{
    int V;    // número de vértices
    list<int> *adj;    // array dinamico de lista de adjacencia

// função que add aresta para o grafo
    void addAresta(int v, int w);

// imprimi cor dos vértices 
    void Colorir();

// função que add aresta
void Grafo::addAresta(int v, int w)

// Atribui cores (a partir de 0) para todos os vértices e imprimi
// A atribuição de cores
void Grafo::Colorir()


#endif //__MDL_H_