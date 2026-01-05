#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#define MAXV 100
#define INF  1e9

int V, E;
double eps;
int MinPts;    

int deg[MAXV];
int to[MAXV][MAXV];
double w[MAXV][MAXV];


int Q[MAXV];          
int inQ[MAXV];        
int density[MAXV];    
int order[MAXV];      
int inOrder[MAXV];    
int orderSize;


int clusterId[MAXV];
int numClusters;


void readGraph(const char* filename) {
    ifstream fin(filename);
    if (!fin.is_open()) {
        cout << "Khong the mo file!" << endl;
        return;
    }

    fin >> V >> E;

    for (int i = 0; i < V; i++)
        deg[i] = 0;

    for (int i = 0; i < E; i++) {
        int u, v;
        double wt;
        fin >> u >> v >> wt;
        
        to[u][deg[u]] = v;
        w[u][deg[u]++] = wt;

        to[v][deg[v]] = u;
        w[v][deg[v]++] = wt;
    }

    fin.close();
}

//--Thuat toan 1- Tim duong di ngan nhat cuc bo
int LSPD(int cp, int inEps[], int isEvent[]){
    double d[MAXV];
    int inQ_LSPD[MAXV];
    int Q_LSPD[MAXV];
    int head = 0, tail = 0;
    
	//Khoi tao d(cp)=0, d() = INF, dday cp vao hang doi 
    for (int i = 0; i < V; i++) {
        d[i] = INF;
        inQ_LSPD[i] = 0;
        inEps[i] = 0;
    }

    d[cp] = 0;
    Q_LSPD[tail++] = cp;
    inQ_LSPD[cp] = 1;

    int dens = 0; //Bien dem so luong phan tu lan can
	//(2) duyet hang doi
    while (head < tail){
    	//(3)
        int p = Q_LSPD[head++];
        inQ_LSPD[p] = 0;
		
		//(4-6)Kiem tra p la su kien-chua co-Them vao N_eps(cp) va tang mat do
        if (isEvent[p] && !inEps[p] && p!=cp ){
            inEps[p] = 1;
            dens++;
        }
		
		//(7)(8)
        if (d[p] < eps){ 
            for (int i = 0; i < deg[p]; i++){
                int q = to[p][i];
                double nd = d[p] + w[p][i];//Tinh khoang cach moi
				
				//Cap nhat khoang cach moi
                if (nd < d[q] && nd <= eps){
                    d[q] = nd;
                    //(11-13)
                    if (!inQ_LSPD[q]){
                        Q_LSPD[tail++] = q;
                        inQ_LSPD[q] = 1;
                    }
                }
            }
        }
    }
    return dens;
}

void pushByDensity(int x, int Q[], int &head, int &tail){
    int i = tail - 1;

    while (i >= head && density[Q[i]] < density[x]) {
        Q[i + 1] = Q[i];
        i--;
    }
    Q[i + 1] = x;
    tail++;
}

//---Thuat toan 2: Tao bang thu tu mat do
void generateDensityOrdering(){
    // 1. Khoi tao lai cac mang
    for (int i = 0; i < V; i++) {
        density[i] = -1;
        inOrder[i] = 0;
        inQ[i] = 0;
    }
    
    orderSize = 0;

    //Tinh nguong mat do
    double threshold = (V > 0) ? log(V) : 0; 
	//(2)
    for (int p = 0; p < V; p++){
    	//(3)
        if (!inOrder[p]){

            int head = 0, tail = 0;
            
            // Chuan bi mang phu tro cho LSPD
            int isEvent[MAXV];
            for (int i = 0; i < V; i++) isEvent[i] = 1;
            int inEpsP[MAXV];

            // Tinh mat do cua p neu co
            if (density[p] == -1) {
                density[p] = LSPD(p, inEpsP, isEvent);
            }

            //Cai Tien: Neu mat do lon hon moi xet
            if (density[p] >= threshold) {
                pushByDensity(p, Q, head, tail); 
                inQ[p] = 1;
            }

            while (head < tail) {
                int q = Q[head++];
                inQ[q] = 0;

                if (density[q] >= threshold){
                    order[orderSize++] = q;
                    inOrder[q] = 1;

                    int inEpsQ[MAXV];
                    LSPD(q, inEpsQ, isEvent);

                    for (int x = 0; x < V; x++){
                        if (inEpsQ[x]) {

                            if (density[x] == -1) {
                                int inEpsX[MAXV];
                                density[x] = LSPD(x, inEpsX, isEvent);
                            }

                            if (!inOrder[x] && !inQ[x] && density[x] >= threshold) {
                                pushByDensity(x, Q, head, tail);
                                inQ[x] = 1;
                            }
                        }
                    }
                }
            }
        }
    }
}

//--Thuat toan 3: Hinh thanh cum
void formClusters() {

    for (int i = 0; i < V; i++)
        clusterId[i] = -1;

    numClusters = 0;

    int isEvent[MAXV];
    for (int i = 0; i < V; i++) isEvent[i] = 1;

    for (int i = 0; i < orderSize; i++) {
        int p = order[i];

        if (clusterId[p] == -1 && density[p] >= MinPts) {
            
            int cid = numClusters++;
            clusterId[p] = cid;

            int Qc[MAXV];
            int head = 0, tail = 0;
            Qc[tail++] = p;

            while (head < tail) {
                int q = Qc[head++];
				
				//tim cac diem bien/ lan can cua q
                int inEpsQ[MAXV];
                LSPD(q, inEpsQ, isEvent);

                for (int x = 0; x < V; x++) {

                    if (inEpsQ[x] && clusterId[x] == -1) {
                        
                        //Mo rong cum: 
                        clusterId[x] = cid;
                        if (density[x] >= MinPts) {
                            Qc[tail++] = x;
                        }
                    }
                }
            }
        }
    }
}

int main() {
    eps = 1.0;
    MinPts = 5;

    readGraph("graph3.txt");
    
    if (V == 0) return 0;

    cout << "--- Thong so ---" << endl;
    cout << "V: " << V << ", E: " << E << endl;
    cout << "Eps: " << eps << ", MinPts: " << MinPts << endl << endl;
	
	int inEpsTest[MAXV];
	int isEvent[MAXV];
	for (int i = 0; i < V; i++) isEvent[i] = 1;
	
	LSPD(0, inEpsTest, isEvent);
	
	cout << "N_eps(P1) = { ";
	for (int i = 0; i < V; i++)
	    if (inEpsTest[i])
	        cout << "P" << i+1 << " ";
	cout << "}\n" << endl;

    generateDensityOrdering();
    cout << "T: ";
    for (int i = 0; i < orderSize; i++)
        cout << "P" << order[i] + 1 << " ";
    cout << endl;

    cout << "Bang thu tu mat do:" << endl;
    for (int i = 0; i < V; i++)
        cout << "P" << i+1 << ": " << density[i] << endl;
    
    
    
    formClusters();
    cout << "\n--- Ket qua hinh thanh cum---" << endl;
    for (int k = 0; k < numClusters; k++) {
        cout << "Cum " << k + 1 << ": { ";
        for (int i = 0; i < V; i++) {
            if (clusterId[i] == k) {
                cout << "P" << i + 1 << " ";
            }
        }
        cout << "}" << endl;
    }

    cout << "Nhieu: { ";
    for (int i = 0; i < V; i++) {
        if (clusterId[i] == -1) {
            cout << "P" << i + 1 << " ";
        }
    }
    cout << "}" << endl;

    return 0;
}
