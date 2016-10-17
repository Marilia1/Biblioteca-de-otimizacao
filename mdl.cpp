<<<<<<< HEAD

=======
#include "mdl.h"
#include <iostream>
#include <cstdlib>


int max(int a, int b)
{
	return a>b? a:b;
}

int min(int a, int b, int c)
{
	return a < b ? (a < c ? a:c): (b < c ? b:c);
}

/*
Seguindo a definição da função:
LCS(Xi,Yj) = 0, se i = 0 ou j = 0;
LCS(Xi,Yj) = LCS(X(i-1), Y(j-1)) + 1, se Xi = Yj;
LCS(Xi,Yj) = max(LCS(Xi, Y(j-1), LCS(X(i-1),Yj)), se Xi != Yj.
*/
int LCS(std::string s1, std::string s2)
{
	int x = s1.length(), y = s2.length();
	int L[x+1][y+1];
	for(int i = 0; i<=x; i++){
		for(int j = 0; j<=y; j++){
			if(!j || !i) L[i][j] = 0;
			else if(s1[i-1] == s2[j-1]) L[i][j] = L[i-1][j-1] + 1;
			else L[i][j] = max(L[i][j-1], L[i-1][j]);
		}
	}
	return L[x][y];
}

/*
Seguindo a definição da função:
SSP(set, n, sum) = false, se sum > 0 e n = 0;
SSP(set, n, sum) = true, se sum = 0;
SSP(set, n, sum) = SSP(set, n-1, sum) || SSP(set, n-1, sum-set[n-1]), caso contrario.
*/
bool SSP(std::vector<int> set, int sum)
{
	int n = set.size();
	bool subset[sum+1][n+1];

	for(int i = 0; i <= n; ++i) subset[0][i] = true;

	for(int i = 0; i <= sum; ++i) subset[i][0] = false;

	for(int i = 1; i <= sum; ++i)
		for(int j = 1; j <= n; ++j){
			subset[i][j] = subset[i][j-1];
			if(i >= set[j-1]) subset[i][j] = subset[i][j] || subset[i - set[j-1]][j-1];
		}

	return subset[sum][n];
}
/*

*/
int DDP(std::string str1, std::string str2)
{
	int m = str1.size();
	int n = str2.size();
    int dp[m+1][n+1];    
    for (int i=0; i<=m; i++)
        for (int j=0; j<=n; j++)      
            if (i==0) dp[i][j] = j;    
            else if (j==0) dp[i][j] = i;            
            else if (str1[i-1] == str2[j-1]) dp[i][j] = dp[i-1][j-1];
			else dp[i][j] = 1 + min(dp[i][j-1], dp[i-1][j], dp[i-1][j-1]); 

    return dp[m][n];
}


MinHeapNode* newNode(char data, unsigned freq)
{
    MinHeapNode* temp = new MinHeapNode();
    temp->left = temp->right = NULL;
    temp->data = data;
    temp->freq = freq;
    return temp;
}
 

MinHeap* createMinHeap(unsigned capacity)
{
    MinHeap* minHeap = new MinHeap();
    minHeap->size = 0;  
    minHeap->capacity = capacity;
    minHeap->array = new MinHeapNode*[minHeap->capacity];
    return minHeap;
}
 

void swapMinHeapNode(MinHeapNode** a, MinHeapNode** b)
{
    MinHeapNode* t = *a;
    *a = *b;
    *b = t;
}
 

void minHeapify(MinHeap* minHeap, int idx)
{
    int smallest = idx;
    int left = 2 * idx + 1;
    int right = 2 * idx + 2;
 
    if (left < minHeap->size &&
        minHeap->array[left]->freq < minHeap->array[smallest]->freq)
      smallest = left;
 
    if (right < minHeap->size &&
        minHeap->array[right]->freq < minHeap->array[smallest]->freq)
      smallest = right;
 
    if (smallest != idx)
    {
        swapMinHeapNode(&minHeap->array[smallest], &minHeap->array[idx]);
        minHeapify(minHeap, smallest);
    }
}
 

int isSizeOne(MinHeap* minHeap)
{
    return (minHeap->size == 1);
}
 

MinHeapNode* extractMin(MinHeap* minHeap)
{
    MinHeapNode* temp = minHeap->array[0];
    minHeap->array[0] = minHeap->array[minHeap->size - 1];
    --minHeap->size;
    minHeapify(minHeap, 0);
    return temp;
}
 

void insertMinHeap(MinHeap* minHeap, MinHeapNode* minHeapNode)
{
    ++minHeap->size;
    int i = minHeap->size - 1;
    while (i && minHeapNode->freq < minHeap->array[(i - 1)/2]->freq)
    {
        minHeap->array[i] = minHeap->array[(i - 1)/2];
        i = (i - 1)/2;
    }
    minHeap->array[i] = minHeapNode;
}
 

void buildMinHeap(MinHeap* minHeap)
{
    int n = minHeap->size - 1;
    int i;
    for (i = (n - 1) / 2; i >= 0; --i)
        minHeapify(minHeap, i);
}
 


int isLeaf(MinHeapNode* root)
{
    return !(root->left) && !(root->right) ;
}
 


MinHeap* createAndBuildMinHeap(char data[], int freq[], int size)
{
    MinHeap* minHeap = createMinHeap(size);
    for (int i = 0; i < size; ++i)
        minHeap->array[i] = newNode(data[i], freq[i]);
    minHeap->size = size;
    buildMinHeap(minHeap);
    return minHeap;
}
 

MinHeapNode* buildHuffmanTree(char data[], int freq[], int size)
{
    MinHeapNode *left, *right, *top;
 
    
    
    MinHeap* minHeap = createAndBuildMinHeap(data, freq, size);
 
    
    while (!isSizeOne(minHeap))
    {
        
        left = extractMin(minHeap);
        right = extractMin(minHeap);
 
        
        
        
        
        top = newNode('$', left->freq + right->freq);
        top->left = left;
        top->right = right;
        insertMinHeap(minHeap, top);
    }
 
    
    return extractMin(minHeap);
}

void printArr(int arr[], int n)
{
    int i;
    for (i = 0; i < n; ++i)
        std::cout << arr[i];
    std::cout << std::endl;
}


void printCodes(MinHeapNode* root, int arr[], int top)
{
    
    if (root->left)
    {
        arr[top] = 0;
        printCodes(root->left, arr, top + 1);
    }
 
    
    if (root->right)
    {
        arr[top] = 1;
        printCodes(root->right, arr, top + 1);
    }
 
    
    
    if (isLeaf(root))
    {
    	std::cout << root->data << ": ";
        printArr(arr, top);
    }
}



void HuffmanCodes(char data[], int freq[], int size)
{
   
   MinHeapNode* root = buildHuffmanTree(data, freq, size);
 
   
   int arr[MAX_TREE_HT], top = 0;
   printCodes(root, arr, top);
}
  

Graph* createGraph(int V, int E)
{
    Graph* graph = new Graph();
    graph->V = V;
    graph->E = E;
 
    graph->edge = new Edge[graph->E];
 
    return graph;
}
 


int find(subset subsets[], int i)
{
    
    if (subsets[i].parent != i)
        subsets[i].parent = find(subsets, subsets[i].parent);
 
    return subsets[i].parent;
}
 


void Union(subset subsets[], int x, int y)
{
    int xroot = find(subsets, x);
    int yroot = find(subsets, y);
 
    
    
    if (subsets[xroot].rank < subsets[yroot].rank)
        subsets[xroot].parent = yroot;
    else if (subsets[xroot].rank > subsets[yroot].rank)
        subsets[yroot].parent = xroot;
 
    
    
    else
    {
        subsets[yroot].parent = xroot;
        subsets[xroot].rank++;
    }
}
 


int myComp(const void* a, const void* b)
{
    Edge* a1 = (Edge*)a;
    Edge* b1 = (Edge*)b;
    return a1->weight > b1->weight;
}
 

void KruskalMST(Graph* graph)
{
    int V = graph->V;
    Edge result[V];  
    int e = 0;  
    int i = 0;  
 
   
    qsort(graph->edge, graph->E, sizeof(graph->edge[0]), myComp);
 
    
    subset *subsets = new subset[V];
 
    
    for (int v = 0; v < V; ++v)
    {
        subsets[v].parent = v;
        subsets[v].rank = 0;
    }
 
    
    while (e < V - 1)
    {
        
        
        Edge next_edge = graph->edge[i++];
 
        int x = find(subsets, next_edge.src);
        int y = find(subsets, next_edge.dest);
 
        
        
        if (x != y)
        {
            result[e++] = next_edge;
            Union(subsets, x, y);
        }
        
    }
 
    
    std::cout << "Following are the edges in the cond MST\n";
    for (i = 0; i < e; ++i)
    	std::cout << result[i].src << " -- " << result[i].dest << " == " << result[i].weight << std::endl;
    return;
}


void shortestPathD(int graph[][V3])
{
	int dist[V3][V3];

	for(int i =0; i<V3 ; ++i)
	{
		for (int j = 0; j < V3; ++j)
		{
			dist[i][j] = graph[i][j];
		}

	}
	for(int i =0; i<V3;++i)
	{
		for(int j =0;j<V3;++j)
		{
			for(int k=0;k<V3;++k)
			{
				if(dist[i][j]!= INF &&
				   dist[j][k]!= INF && 
					dist[i][j]+dist[j][k]<dist[i][k])
				{
					dist[i][k] = dist[i][j] + dist[j][k];
				}
			}
		}
	}
	  for (int i = 0; i < V3; i++)
    	{
		  for (int j = 0; j < V3; j++)
	        {
	            if (dist[i][j] == INF)
	            {
	                 std::cout<<" INF" ;
	            }
	            else
	            {
	                std::cout<<" " << dist[i][j];
	            }
	        }
        std::cout<<std::endl;
    		}
		}

Box::Box(){
	h =0;
	w = 0;
	d = 0;
	calcBase();
}
Box::Box(int i,int j,int k){

	h = i;
	w = j;
	d = k;
	calcBase();
}

int Box::calcBase(){
	base = w * d;
	return base;
}

int min1 (int x, int y)
{ return (x < y)? x : y; }

int max1 (int x, int y)
{ return (x > y)? x : y; }

int maxHeight(std::vector<Box> &x ,int n){
	std::vector<Box> rot(3*n);
	int index = 0;
	for(int i = 0; i<n ;i++){
		rot[index]=x[i];
		index++;

	  	rot[index].h = x[i].w;
     	rot[index].d = max1(x[i].h, x[i].d);
      	rot[index].w = min1(x[i].h, x[i].d);
      	rot[index].calcBase();
      	index++;
 
    
      	rot[index].h = x[i].d;
      	rot[index].d = max1(x[i].h, x[i].w);
      	rot[index].w = min1(x[i].h, x[i].w);
      	rot[index].calcBase();
      	index++;	
	}

	n = 3*n;
	std::sort(rot.begin(),rot.end(),EntityComp());
	
	int msh[n];
   for (int i = 0; i < n; i++ )
      msh[i] = rot[i].h;
 
   
   for (int i = 1; i < n; i++ )
      for (int j = 0; j < i; j++ )
         if ( rot[i].w < rot[j].w &&
              rot[i].d < rot[j].d &&
              msh[i] < msh[j] + rot[i].h
            )
         {
              msh[i] = msh[j] + rot[i].h;
         }
 
 
 
   int max = -1;
   for ( int i = 0; i < n; i++ )
      if ( max < msh[i] )
         max = msh[i];
 
   return max;
}


int dictionaryContains(std::string word)
{
    std::string dictionary[] = {"mobile","samsung","sam","sung","man","mango",
                           "icecream","and","go","i","like","ice","cream"};
    int size = sizeof(dictionary)/sizeof(dictionary[0]);
    for (int i = 0; i < size; i++)
        if (dictionary[i].compare(word) == 0)
           return true;
    return false;
}

std::string lastWord(std::string s)
{
    std::string x = " ";
    std::string y;

    for(int i = s.size()-2; i>0;--i)
    {
        if(s[i] == x[0] || i==0 )
        {
       
            y= s.substr(i,s.size()-i);
            break;
        }
    }
    return y;
}

bool wordBreak(std::string str)
{
     
     std::string b;
     std::string a;
     std::string c;   
    int size = str.size();
    if (size == 0)   return true;
 
    
    bool wb[size+1];
    memset(wb, 0, sizeof(wb)); 
 
    for (int i=1; i<=size; i++)
    {
        
        if (wb[i] == false && dictionaryContains( str.substr(0, i) ))
        {
            wb[i] = true;
            b.append(" ");
            b.append(str.substr(0,i));
            b.append(" ");
         
        }
        
        if (wb[i] == true)
        {
          
         
 
            for (int j = i+1; j <= size; j++)
            {
                
                if (wb[j] == false && dictionaryContains( str.substr(i, j-i) ))
                {
                  
                    wb[j] = true;
                    c =str.substr(i,j-i);
                    a = lastWord(b);
                    a.erase(0,1);
                    if(a.size()>1){
                    a.pop_back();
                  }
                  

                    if(c.find(a)!= std::string::npos)
                    {
                        if(!c.compare(a))
                        {
                           
                          b.append(str.substr(i,j-i));
                          b.append(" ");       
                        }
                        else 
                        {
                           
                         b.erase(b.size()-(a.size()+1),a.size());
                         b.append(str.substr(i,j-i));
                         b.append(" ");    
                     }
                    }
                    if(c.find(a)==std::string::npos)
                    {
                        
                    b.append(str.substr(i,j-i));
                    b.append(" ");    
                    }
                   
                }
            }
 
                
                
            }
        }
    
 
 
    std::cout<< b <<std::endl;
 
   
    return false;
}

int minDistance(int dist[], bool sptSet[])
{
   // Initialize min value
   int min = INT_MAX, min_index;
  
   for (int v = 0; v < V1; v++)
     if (sptSet[v] == false && dist[v] <= min)
         min = dist[v], min_index = v;
  
   return min_index;
}

int printSolution(int dist[], int n)
{
   std::cout<<"Vertex   Distance from Source"<<std::endl;
   for (int i = 0; i < V1; i++)
      std::cout<< i <<" " <<dist[i]<<std::endl;;
}

void dijkstra(int graph[V1][V1], int src)
{
     int dist[V1];     
  
     bool sptSet[V1]; 
  
    
     for (int i = 0; i < V1; i++)
        dist[i] = INT_MAX, sptSet[i] = false;
  
    
     dist[src] = 0;
  
     
     for (int count = 0; count < V1-1; count++)
     {
       
       int u = minDistance(dist, sptSet);
  
       
       sptSet[u] = true;
  
     
       for (int v = 0; v < V1; v++)
  
         
         if (!sptSet[v] && graph[u][v] && dist[u] != INT_MAX 
                                       && dist[u]+graph[u][v] < dist[v])
            dist[v] = dist[u] + graph[u][v];
     }
  
    
     printSolution(dist, V1);
}

int minKey(int key[], bool mstSet[])
{
  
   int min = INT_MAX, min_index;
 
   for (int v = 0; v < V2; v++)
     if (mstSet[v] == false && key[v] < min)
         min = key[v], min_index = v;
 
   return min_index;
}

int printMST(int parent[], int n, int graph[V2][V2])
{
   std::cout<<"Edge   Weight"<<std::endl;
   for (int i = 1; i < V2; i++)
    std::cout<<parent[i]<<" " << i << " " <<graph[i][parent[i]]<<std::endl;
}

void primMST(int graph[V2][V2])
{
     int parent[V2]; 
     int key[V2];  
     bool mstSet[V2];  
 
    
     for (int i = 0; i < V2; i++)
        key[i] = INT_MAX, mstSet[i] = false;
 
     
     key[0] = 0;     
     parent[0] = -1; 
 
   
     for (int count = 0; count < V2-1; count++)
     {
        
        int u = minKey(key, mstSet);
 
     
        mstSet[u] = true;
 
       
        for (int v = 0; v < V2; v++)
 
          
          if (graph[u][v] && mstSet[v] == false && graph[u][v] <  key[v])
             parent[v]  = u, key[v] = graph[u][v];
     }
 
    
     printMST(parent, V2, graph);
}


//-----------------Programação Dinamica------------------
//Os algoritmos trabalhados na programação dinamica sera:
//Triângulo de Pascal
//Troco em moedas
//Mochila binária

/*O primeiro programa é o triângulo de Pascal onde criamos uma função 
que recebem como parametro as combinações que desejamos procurar no
triangulo, com ela sabaremos qual numero que pertence aquela linha e coluna */



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

>>>>>>> refs/remotes/origin/pr/2
