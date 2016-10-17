<<<<<<< HEAD
#include <iostream>
#include <list>
#include <algorithm>
#include <limits.h>
#include <cstdio>
#include <vector>
#include <string>


=======
<<<<<<< HEAD
=======
<<<<<<< HEAD
>>>>>>> refs/remotes/origin/master
>>>>>>> refs/remotes/origin/master
#ifndef __MDL_H__
#define __MDL_H__
#define INF INT_MAX
#define V1 9
#define V2 5
#define V3 4
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




//--------------- Pogramação dinamica ------------------
/*
Retorna o maior entre dois inteiros.
Parametros: Numeros a serem comparados.
*/
int max(int, int);

int min(int, int, int);

/*
Aluno: Lucas Ferreira Lima - 11406537
Encontra a maior subsequência comum a duas sequências X e Y, i.e uma 
subsequência de uma sequência X de caracteres pode ser obtida removendo alguns 
caracteres dessa sequência.

Parametros: Subsequêcias a serem comparadas.
Retorno: Valor da maior subsequencia.
*/
int LCS(std::string, std::string);

/*
Aluno: Lucas Ferreira Lima - 11406537
Encotra em um dado conjunto de inteiros positivos se existe algum subconjunto 
cujo a soma é sum.

Parametros: Vector com conjunto de inteiros positivos e soma que deseja ser
			comparada.
Retorno: Numero de subconjuntos encontrados.
*/
bool SSP(std::vector<int> , int);

/*
Aluno: Lucas Ferreira Lima - 11406537
Encontra o número mínimo de operações necessárias para transformar um string no outro.
Sendo essas oeprações: adição, remoção e substituíção

Parametros: Strings a serem comparadas.
Retorno: Numero de operções nescessárias.
*/
int DDP(std::string, std::string);

//--------------- Algoritmos gulosos ------------------
/*
Aluno: Lucas Ferreira Lima - 11406537
Huffman Coding
*/

/*
Nó da MinHeap. 
Contém um campo para o dado (um dos caracteres de input),
frequênica do dado,
e aponta para os filhos da direita e esquerta.
*/
typedef struct MinHeapNode
{
    char data;
    unsigned freq;
    MinHeapNode *left, *right;
} MinHeapNode;

/*
Uma MinHeap. 
Contém o tamanho da MinHeap,
capacidade da mesma,
e um array de todos os filhos.
*/
typedef struct MinHeap
{
    unsigned size;
    unsigned capacity;
    MinHeapNode **array;
} MinHeap;

/*
Função que aloca um nova minHeap.
Parametro: Um caractere e sua frequencia.
Retorno: MinHeap com dado caractere e frequencia.
*/
MinHeapNode* newNode(char data, unsigned freq);

/*
Função que cria uma MinHeap de dada capacidade.
Parametro: Capacidade da MinHeap a ser craida.
Retorno: MinHeap com dada capacidade.
*/
MinHeap* createMinHeap(unsigned capacity);

/*
Troca dois nó de uma MinHeap.
Parametro: Dois nó de uma MinHeap.
Retorno: Não há.
*/
void swapMinHeapNode(MinHeapNode** a, MinHeapNode** b);
/*
Função padrão de minHeapify
Parametro: MinHeap e index.
Retorno: Não há
*/
void minHeapify(MinHeap* minHeap, int idx);
/*
Checa se o tamanho da heap é 1 ou não.
Paramentro: Uma MinHeap.
Retorno: 0 se o tamnho for diferente de 1, diferente de 0 se o
	tamanho for 1.
*/
int isSizeOne(MinHeap* minHeap);
/*
Extrai o menor nó de menor valor da MinHeap.
Parametro: Uma MinHeap.
Retorno: O nó de menor valor da MinHeap.
*/
MinHeapNode* extractMin(MinHeap* minHeap);

/*
Insere um novo nó à uma MinHeap.
Parametro: Um nó a ser adcionado na MinHeap e respectiva MinHeap.
Retorno: Não há.
*/
void insertMinHeap(MinHeap* minHeap, MinHeapNode* minHeapNode);

/*
Função para construir MinHeap.
Parametro: MinHeap á qual será construida.
Retorno: Não há.
*/
void buildMinHeap(MinHeap* minHeap);

/*
Checa se o nó é folha
Paramentro: Uma MinHeap.
Retorno: 0 se o nó não for folha, diferente de 0 se o nó for folha.
*/
int isLeaf(MinHeapNode* root);

/*
Cria uma MinHeap de capacidade igual a o tamanho e insere todos os caracteres
do array data na MinHeap. Inicializa o tamanha da MinHeap com a capacidade.
Parametro: array de dados que serão inseridos na MinHeap, frequêcia dos dados,
e o tamanho.
Retorno: MinHeap criada dito as especificações da função.
*/
MinHeap* createAndBuildMinHeap(char data[], int freq[], int size);

/*
Constrói a árvore de Huffman 
Parametro: array de dados, frequencia dos dados e o tamanho.
Retorno: Huffman coding
*/
MinHeapNode* buildHuffmanTree(char data[], int freq[], int size);

/*
Função principal que constrói e imprime Huffman code.
Parametro: Array de dados, frequencia dos dados e tamanho.
Retorno: Não há.
*/ 
void HuffmanCodes(char data[], int freq[], int size);

/*
Imprime na tela o codigo de Huffman da raíz da Àrvore de Huffman. Usa o array
arr para aramazenar os códigos.
Parametro: MinHeap que é a Árvore de Huffman, array para o armazenamento dos
	códigos e o topo.
Retorno: Não há.
*/
void printCodes(MinHeapNode* root, int arr[], int top);

/*
Imprime uma array de tamnho n nao tela.
Parametro: array e tamanho.
Retorno: Não há.
*/
void printArr(int arr[], int n);
//----------------------------------------------------------------------------------------------
/*
Aluno: Lucas Ferreira Lima - 11406537
Kruskal MST
*/

/*
Estrutura para representar aresta de um grafos com peso.
Contém a fonte, o destino e o peso.
*/
// a structure to represent a weighted edge in graph
typedef struct Edge
{
    int src, dest, weight;
} Edge;
/*
Estrutura para representar um grafo conectado, indireto e com peso.
Contém numero de vertices (V), numero de arestas (E), e um array de arestas.
*/ 
typedef struct Graph
{
    int V, E;
    Edge* edge;
} Graph;

typedef struct subset
{
    int parent;
    int rank;
} subset;

/*
Cria grafo com V vertices e E arestas.
Parametro: Numero de vertices e de arestas.
Retorno: Grafo descrito na função.
*/
Graph* createGraph(int V, int E);

/*
Faz a união de dois dois sets. (Usa união por rank)
Parametro: Um subset, e dois set.
Retorno: Não há.
*/
void Union(subset subsets[], int x, int y);

/*
Compara duas arestas de acordo com seu peso. Usa qsort com algoritimo de 
ordenação.
Parametro: Duas arestas a serem comparada.
Retorno: 
*/
int myComp(const void* a, const void* b);

/*
Função principal para construir MST usando algoritmo de Kruskal.
Parametro: Grafo que se deseja construir a MST.
Retorno: Não há.
*/
void KruskalMST(Graph* graph);

/*Função que recebe como parametro o grafo 
como um array multidimensional de inteiros
e imprime uma "tabela" com as menores distancias*/
void shortestPathD(int graph[][V3]);	

	
/*classe que representa uma caixa*/
class Box{
public:
	int h,w,d;
	int base;
	

Box();
Box(int i,int j,int k);

int calcBase();
	

};


/*Struct para ser usada na função de sort*/
struct EntityComp {
 
  bool operator()(const Box& s1, const Box& s2) const {
      
          return s1.base > s2.base;
      }
};

/*Função para retornar o minimo entre dois numeros,
parametros:inteiros*/
int min1 (int x, int y);

/*Função para retornar o max entre dois numeros,
parametros: inteiros*/
int max1 (int x, int y);
	
/*Função para retornar a maior altura que podemos
empilhar as caixas, parametros: vector de caixas 
e numero de caixas*/
int maxHeight(std::vector<Box> &x ,int n);

/*Função para checar se uma palavra está no dicionário ou não,
entrada uma string retorna um inteiro 1 se estiver e 
zero se nao estiver*/
int dictionaryContains(std::string word);

/*função que recebe como entrada uma string
e retorna a ultima palavra dela
(considerando que a string comece e termine com 
espaço em branco)*/
std::string lastWord(std::string s);

/*Função que retorna false se nao
der pra dividir(de acordo com o dicionário)
 imprime uma opção de repartição*/
bool wordBreak(std::string str);
	
/*Função para achar o vertice com a menor distancia
do conjunto de vertices que nao estao na arvore do menor caminho
entrada : um array de inteiros com as distancias e um array de booleans
dziendo se está contido ou nao*/
int minDistance(int dist[], bool sptSet[]);

/*Função para imprimir a "tabela" de distancia
entrada: um array de inteiros com as distancias e
o numero de vertices */
int printSolution(int dist[], int n);

/*Função que acha o caminho mais curto usando o algoritmo
de dijkstra para um grafo representado por matriz de adjacencia 
entrada o grafo como um array multidimensional de inteiros e a 
origem como um inteiro*/
void dijkstra(int graph[V1][V1], int src);

/*Função para achar o vertice com a menor 
chave do conjunto de vertices não incluídos na Arvore 
entrada um array de inteiros com as chaves e um array de 
booleans com true se estiver contido e false se nao
estiver contido */
int minKey(int key[], bool mstSet[]);

/*Função para imprimir a arvore armazenada em parent[]
entrada a arvore como um array, o numero de vertices como um inteiros e
o grafo como um array de inteiros multidimensional*/
int printMST(int parent[], int n, int graph[V2][V2]);

/*Função para construir e imprimir o grafo
usando a matriz de adjacencia, entrada: o grafo 
como um array multidimensional de inteiros*/
void primMST(int graph[V2][V2]);

	
	
#endif //__MDL_H_
>>>>>>> refs/remotes/origin/pr/2
>>>>>>> refs/remotes/origin/master
