//-----------------Programação Dinamica------------------
//Os algoritmos trabalhados na programação dinamica sera:
//Triângulo de Pascal
//Troco em moedas
//Mochila binária

/*O primeiro programa é o triângulo de Pascal onde criamos uma função 
que recebem como parametro as combinações que desejamos procurar no
triangulo, com ela sabaremos qual numero que pertence aquela linha e coluna */

#include <iostream>

int min(int a, int b);

int Coeficiente(int lin, int col)
{
    int C[lin+1][col+1];
    int i, j;

    for (i = 0; i <= lin; i++)
    {
        for (j = 0; j <= min(i, col); j++)
        {
            // caso base
            if (j == 0 || j == i)
                C[i][j] = 1;

            // Calcular o valor usando os valores armazenados anteriormente
            else
                C[i][j] = C[i-1][j-1] + C[i-1][j];
        }
    }

    return C[lin][col];
}

// retorna min de dois inteiros
int min(int a, int b)
{
    return (a<b)? a: b;
}
//-----------------Fim do primeiro programa-------------

/*O segundo programa é o troco em moedas, esse programa nos diz quantas 
possibilidades de combinação para um certo valos de troco podemos ter 
de acordo com os valores de moedas disponíveis */


int count( int S[], int m, int n )
{
    int i, j, x, y;

    // Precisaremos de n+1 linhas para contruir uma tabela
    // caso base 0 (n = 0)
    int tabela[n+1][m];

   //Caso o n = 0 preencha as entrada para 0
    for (i=0; i<m; i++)
        tabela[0][i] = 1;

    // Preenche o resto da tabela com as entradas, de baixo para cima
    for (i = 1; i < n+1; i++)
    {
        for (j = 0; j < m; j++)
        {
            // Contagem de soluções incluindo o S[j]
            x = (i-S[j] >= 0)? tabela[i - S[j]][j]: 0;

            // Contagem de soluções excluindo o S[j]
            y = (j >= 1)? tabela[i][j-1]: 0;

            // Total de contagens
            tabela[i][j] = x + y;
        }
    }
    return tabela[n][m-1];
}
//------------------Fim do segundo programa-----------------

/*O terceiro programa é o de mochila binária, onde a ideia é que temos uma mchila que 
tem a capacidade para suportar um peso n, e temos objetos com diferentes pesos, e queremos
saber quais itens a mochila aguenta*/

int max(int a, int b) { return (a > b)? a : b; }

// Retrona o valor máximo que pode ser colocado na mochila com uma capacidade C
int knapSack(int C, int pitem[], int val[], int n)
{
   int i, peso;
   int tabela[n+1][C+1];

   // Construindo uma tabela[][] de baixo pra cima
   for (i = 0; i <= n; i++)
   {
       for (peso = 0; peso <= C; peso++)
       {
           if (i==0 || peso==0)
               tabela[i][peso] = 0;
           else if (pitem[i-1] <= peso)
                 tabela[i][peso] = max(val[i-1] + tabela[i-1][peso-pitem[i-1]],  tabela[i-1][peso]);
           else
                 tabela[i][peso] = tabela[i-1][peso];
       }
   }

   return tabela[n][C];
}
//------------------------Fim do terceiro programa--------------

//--------------------Algoritmos gulosos----------------
//Os algoritmos trabalhados nos gulosos serão:
//Mochila fracionária
//Coloração de grafos

/*O primeiro programa é a mochila fracionária, é basicamente a mesma ideia do algoritmo 
mostrado acima, porem dessa vez caso um item tenha um peso que a mochila não suporta
apenas parte dele será add obtendo o valor total da capacidade C da mochila*/

struct Item
{
    int valor, peso;
 
    // Construtor
    Item(int valor, int peso) : valor(valor), peso(peso)
    {}
};
 
// Função que compara os intes de acordo com a relação do seu valor pelo peso
bool cmp(struct Item a, struct Item b)
{
    double r1 = (double)a.valor / a.peso;
    double r2 = (double)b.valor / b.peso;
    return r1 > r2;
}
 
// principal função que resolve o problema
double fractionalKnapsack(int C, struct Item arr[], int n)
{
    //Ordenar itensbaseado na relação
    sort(arr, arr + n, cmp);
    int pesoatual = 0;  // peso atual na mochila
    double valorfinal = 0.0; // resultado (valor na mochila)
 
    // percorrendo todos os itens
    for (int i = 0; i < n; i++)
    {
        // caso o item não ultrapasse a capacidade da mochila, add
        if (pesoatual + arr[i].peso <= C)
        {
            pesoatual += arr[i].peso;
            valorfinal += arr[i].valor;
        }
 
        // Caso não possamos add o item todo, add parte dele
        else
        {
            int resta = C - pesoatual;
            valorfinal += arr[i].valor * ((double) resta / arr[i].peso);
            break;
        }
    }
 
    // retorna o valor final
    return valorfinal;
}
//---------------------Fim do primeiro programa---------------------

/*O segundo programa é Coloração de grafos,onde o problema é, dado m cores, 
temos que encontrar uma maneira de colorir os vértices de um grafo de modo que 
não tenha dois vértices adjacentes com a mesma cor. Esse algoritmo não garante 
utilizar cores mínimas, mas garante um limite superior do número de cores. 
O algoritmo básico nunca usa mais do que d cores + 1, onde d é o grau máximo 
de um vértice no gráfico dado. */

#include <list>

//Classe que representa um grafo não direcionavel
class Grafo
{
    int V;    // número de vértices
    list<int> *adj;    // array dinamico de lista de adjacencia
public:
    // Construtor e desconstrutor
    Grafo(int V)   { this->V = V; adj = new list<int>[V]; }
    ~Grafo()       { delete [] adj; }

    // função que add aresta para o grafo
    void addAresta(int v, int w);

    // imprimi cor dos vértices 
    void Colorir();
};

void Grafo::addAresta(int v, int w)
{
    adj[v].push_back(w);
    adj[w].push_back(v);  
}

// Atribui cores (a partir de 0) para todos os vértices e imprimi
// A atribuição de cores
void Grafo::Colorir()
{
    int resultado[V];

    // Atribui a primeira cor para o primeiro vértice
    resultado[0]  = 0;

    // Inicializar restantes V-1 vértices como não atribuído
    for (int u = 1; u < V; u++)
        resultado[u] = -1;  // Não atribui cor a u

    // Array temporario para add as cores. True
     // Valor disponível [cr] significaria que a cor cr é
     // Atribuído a um de seus vértices adjacentes
    bool disponivel[V];
    for (int cr = 0; cr < V; cr++)
        disponivel[cr] = false;

    // Atribuir cores para restantes V-1 vértices
    for (int u = 1; u < V; u++)
    {
        //processa todos os vértices adjacentes e coloca suas cores como indisponivel
        list<int>::iterator i;
        for (i = adj[u].begin(); i != adj[u].end(); ++i)
            if (resultado[*i] != -1)
                disponivel[resultado[*i]] = true;

        // Encontra a primeira cor disponivel
        int cr;
        for (cr = 0; cr < V; cr++)
            if (disponivel[cr] == false)
                break;

        resultado[u] = cr; // Atribui cor encontrada

        // Volta os valores para falso na proxima interação
        for (i = adj[u].begin(); i != adj[u].end(); ++i)
            if (resultado[*i] != -1)
                disponivel[resultado[*i]] = false;
    }

    // 
    for (int u = 0; u < V; u++)
        cout << "Vertice " << u << " --->  Cor "
             << resultado[u] << endl;
}
//------------------Fim do segundo programa-----------------