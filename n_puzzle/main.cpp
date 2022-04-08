#include <bits/stdc++.h>
#include <unistd.h>
#include <chrono>
#define N 4
using namespace std;
using namespace std::chrono;

int dir_x[4] = { 1, 0, -1, 0 };
int dir_y[4] = { 0, -1, 0, 1 };

class Node {
    public:
        bool game_over();
        void fill_random();
        void fill_shuffle(int s);
        void print();
        void calculate_distance();
        void calculate_misplaced();
        string hash();

        int find_X_position();
        bool is_solvable();

        Node(int b[N][N]);

        int board[N][N];
        int empty_x = N-1;
        int empty_y = N-1;
        int heuristic = 0;
        int moves = 0;
};

Node::Node(int b[N][N]){
    for (int i=0; i< N;i++){
		for(int j=0;j< N;j++){
			board[i][j] = b[i][j];
		}
	} 
}

string Node::hash(){
    string s = "";
    for (int i=0; i< N;i++){
        for(int j=0;j< N;j++){
            s += to_string(board[i][j]);
        }
    }
    return s;
}

bool Node::game_over(){
    for (int i=0; i< N;i++){
        for(int j=0;j< N;j++){
            if (i == N-1 && j == N-1){
                continue;
            }
            if (board[i][j] != i * N + j + 1){
                return false;
            }
        }
    }
    return true;
}

void Node::calculate_distance(){
    heuristic = 0;
    for (int i=0; i< N;i++){
        for(int j=0;j< N;j++){
            int val = board[i][j];
            if (val == 0){
                continue;
            }
            int x = val/N;
            if (val % N == 0){
                x = x - 1;
            }
            int y = val - x * N - 1;
            heuristic += abs(x - i) + abs(y - j);
        }
    }
}

void Node::calculate_misplaced(){
    heuristic = 0;
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++)
        if (board[i][j] && board[i][j] != i * N + j + 1)
           heuristic++;
}

void Node::fill_shuffle(int shuffles){
    for (int i=0; i< N; i++){
        for(int j=0; j< N; j++){
            if (i == N-1 && j == N-1){
                continue;
            }
            board[i][j] = i * N + j + 1;
        }
    }
    empty_x = empty_y = N-1;
    for (int i = 0; i<shuffles; i++){
        int rand_ind = rand() % 4;
        int x = empty_x + dir_x[rand_ind];
        int y = empty_y + dir_y[rand_ind];

        if (x < 0 || x >= N){
            continue;
        }
        if (y < 0 || y >= N){
            continue;
        }

        board[empty_x][empty_y] = board[x][y];
        board[x][y] = 0;

        empty_x = x;
        empty_y = y;
    }
}

void Node::fill_random(){
    vector<int> cells;
    for (int i = 1; i < N*N; i++){
        cells.push_back(i);
    }
    auto rng = std::default_random_engine {};
    std::shuffle(std::begin(cells), std::end(cells), rng);
    for (int i=0; i< N;i++){
        for(int j=0;j< N;j++){
            if (i == N-1 && j == N-1){
                continue;
            }
            int rand_ind = rand() % cells.size();
            board[i][j] = cells[rand_ind];
            cells.erase(cells.begin() + rand_ind);
        }
    }
}

void Node::print(){
    int r;
    // cout<<heuristic<< " "<<moves<<endl;
    for(int i = 0; i < N; i++) {
        cout << "+----+----+----+----+\n";
        for( int j = 0; j < N; j++) {
            r = board[i][j];
            cout << "| ";
            if(r < 10) cout << " ";
            if(!r) cout << "  ";
            else cout << r << " ";
        }
        cout << "|\n";
    }
    cout << "+----+----+----+----+\n";
}

struct compare{
    public:
    bool operator()(Node a,Node b) {
        if (a.heuristic + a.moves == b.heuristic + b.moves){
            return a.moves < b.moves;
        }
        return a.heuristic + a.moves > b.heuristic + b.moves; 
    }
};

int get_inv_count(int arr[]){
    int inv_count = 0;
    for (int i = 0; i < N * N - 1; i++){
        for (int j = i + 1; j < N * N; j++){
            if (arr[j] && arr[i] && arr[i] > arr[j])
                inv_count++;
        }
    }
    return inv_count;
}
 
int Node::find_X_position(){
    for (int i = N - 1; i >= 0; i--){
        for (int j = N - 1; j >= 0; j--){
            if (board[i][j] == 0){
                return N - i;
            }
        }
    }
    return 0;
}

bool Node::is_solvable(){
    int invCount = get_inv_count((int*)board);
    if (N & 1)
        return !(invCount & 1);
    else{
        int pos = find_X_position();
        if (pos & 1)
            return !(invCount & 1);
        else
            return invCount & 1;
    }
}

int main(){
    // int board[N][N] = {{1, 5, 2},{4, 3, 0},{7, 8, 6}};
    int board[N][N] = {0};
    Node root(board);
    // root.empty_x = 1;
    // root.empty_y = 2;
    root.fill_shuffle(500);
    root.calculate_distance();
    if (root.is_solvable()){
        root.print();
        cout<<"*****************Solvable*****************\n";
        sleep(1);
    }

    map<string, int> mp;
    mp.insert({root.hash(), 1});

    priority_queue<Node,vector<Node>,compare> q;
    q.push(root);

    auto start = high_resolution_clock::now();

    while (q.size() > 0){
        Node curr = q.top();
        q.pop();

        if (curr.game_over() == true){
            curr.print();
            break;
        }

        for (int d = 0; d <4; d++){
            int x = curr.empty_x + dir_x[d];
            int y = curr.empty_y + dir_y[d];

            if (x < 0 || x >= N){
                continue;
            }
            if (y < 0 || y >= N){
                continue;
            }

            Node next(curr.board);
            next.board[curr.empty_x][curr.empty_y] = curr.board[x][y];
            next.board[x][y] = 0;
            next.empty_x = x;
            next.empty_y = y;
            next.moves = curr.moves +1;
            next.calculate_distance();
            if (mp[next.hash()] == 0){
                mp[next.hash()] = 1;
                q.push(next);
            }
        }
    }

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    cout<<"Time taken(s): "<<duration.count()<<endl;
    return 0;
}
