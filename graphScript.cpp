#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <fstream>
#include <set>
#include <algorithm>

using namespace std;

void printVec(vector<int> & vec, ostream &fout);
bool readPref(vector<int> & p, int & n, ifstream &fin);
void preamble(ofstream &fout);
void printPGraph(vector<int> & p, ofstream &fout);
void printQGraph(vector<int> & p, ofstream &fout);
void printDoubleVec(vector<vector<int> > vec);
void printVertexEdges(vector<vector<pair<int, int> > > s);
bool pairEq(pair<int, int> a, pair<int, int> b);
bool pairComp(pair<int, int> a, pair<int, int> b);
int coloring(vector<vector<int> > &adj, vector<int> &colors);
int greedyColoring(vector<vector<int> > &adj, vector<int> &colors, vector<int> &perm, int best);
pair<vector<vector<int> > , vector<pair<int, int> > > printQColoring(vector<int> & p, int &red);
int index(pair<int, int> p, vector<pair<int, int> > &positive);
vector<vector<int> > indetString(vector<int> p, vector<int> &colors, vector<pair<int, int> > &positive, int &k);
void printIndetString(vector<vector<int> > &w);

int main()
{
    int n = 0;
    vector<int> p;
    ifstream fin;
    fin.open("pref.txt");

    ofstream fout;
    fout.open("scriptTest.tex");
    preamble(fout);

    int max = 0;
    int k = 0;
    int red = 0;
    while(readPref(p, n, fin))
    {
        //printPGraph(p, fout);
        printQGraph(p, fout);

        pair<vector<vector<int> > , vector<pair<int, int> > > temp = printQColoring(p, red);
        //printDoubleVec(temp.first);
        vector<int> colors((int)temp.first.size(), 0);
        k = coloring(temp.first, colors);
        indetString(p, colors, temp.second, k);
        p.push_back(k);
        printPGraph(p, fout);
        if(k > red)
        {
            printVec(p, cout);
            cout << "Num red edges: " << red << endl;
            cout << "Chromatic number " << k << endl;
        }

        if(k > max)
        {
            max = k;
        }
    }
    cout << "The largest chromatic number is: " << max << endl;

    fout << "\\end{document}" << endl;
    return 0;
}

void printVec(vector<int> & vec, ostream &fout)
{
    for(int i = 0; i < (int)vec.size(); i++)
    {
        fout << vec[i] << " ";
    }
}

bool readPref(vector<int> & p, int & n, ifstream &fin)
{
    fin >> n;
    int e = 0;
    p.clear();
    p.push_back(n);
    for(int i = 1; i < n; i++)
    {
        if(!fin)
        {
            return false;
        }
        fin >> e;
        p.push_back(e);
    }

    //cout << "Prefix array: " << endl;
    //printVec(p, cout);
    //cout << endl;
    return true;
}

void preamble(ofstream &fout)
{
    fout << "\\documentclass{article}" << endl
    << "\\usepackage[dvipsnames]{xcolor}" << endl
    << "\\usepackage{tikz}" << endl << endl
    << "\\tikzstyle{vertex}=[circle, thick, draw=black, inner sep = 2pt]" << endl
    << "\\tikzstyle{positive}=[ultra thick, color=blue]" << endl
    << "\\tikzstyle{negative}=[ultra thick, color=red]" << endl << endl
    << "\\begin{document}" << endl << endl;
}

void printPGraph(vector<int> & p, ofstream &fout)
{
    int k = p.back();
    p.pop_back();
    //fout << "\\begin{figure}" << endl
    //<< "\\centering" << endl
    fout << "\\begin{tikzpicture}[scale=.5]" << endl;

    int n = (int)p.size();
    fout << "\\foreach \\i in {0,1,2,...," << n-1 << "}" << endl
    << "\\node[vertex] (\\i) at ({90 + " << (double)360/n << "*(\\i)}:3){\\i};" << endl << endl;

    fout << "\\node at (270:4) {k is " << k << ". $\\mathcal{P}$ for $y = ";
    printVec(p, fout);
    fout << "$};" << endl;

    for(int i = 1; i < n; i++)
    {
        for(int j = 1; j <= p[i]; j++)
        {
            if(j + i <= n)
            fout << "\\draw[positive] (" << j - 1 << ")--("
            << j + i - 1 << ");" << endl;
        }
        if(p[i] + i + 1 <= n)
        {
            fout << "\\draw[negative] (" << p[i] << ")--("
            << p[i] + i << ");" << endl;
        }

    }

    fout << "\\end{tikzpicture}" << endl;
    //<< "\\caption{$\\mathcal{P}$ graph for ";
    //printVec(p, fout);
    //fout << "}" << endl
    //<< "\\end{figure}" << endl << endl;
}

void printQGraph(vector<int> & p, ofstream &fout)
{
    int n = (int)p.size();

    vector<int> temp(n, 0);
    vector<vector<int> > aPlus(n, temp);
    vector<vector<int> > aMinus = aPlus;

    for(int i = 1; i < n; i++)
    {
        for(int j = 1; j <= p[i]; j++)
        {
            if(j + i <= n)
            {
                aPlus[j - 1][j + i - 1] = 1;
                aPlus[j + i - 1][j - 1] = 1;
            }
        }
        if(p[i] + i + 1 <= n)
        {
            aMinus[p[i]][p[i] + i] = 1;
            aMinus[p[i] + i][p[i]] = 1;
        }
    }
    //printDoubleVec(aPlus);
    //printDoubleVec(aMinus);

    vector<pair<int, int> > stuff;
    vector<vector<pair<int, int> > > s(n, stuff);
    pair <int, int> edge;

    for(int i = 0; i < n; i++) //for each vertex in P
    {
        for(int j = 0; j < n; j++)
        {
            if(aMinus[i][j]) //if it's connected by a red edge
            {
                for(int k = 0; k < n; k++)
                {
                    if(aPlus[j][k])
                    {
                        if(j < k)
                        {
                            edge = make_pair(j, k);
                        }
                        else
                        {
                            edge = make_pair(k, j);
                        }
                        s[i].push_back(edge);
                    }
                }
            }
        }
        sort(s[i].begin(), s[i].end(), pairComp);
        vector<pair<int, int> >::iterator it;
        it = unique(s[i].begin(), s[i].end(), pairEq);
        s[i].resize(it - s[i].begin());
    }
    //printVertexEdges(s);

    vector<pair<int, int> > positive;
    for(int i = 0; i < n; i++)
    {
        for(int j = i; j < n; j++)
        {
            if(aPlus[i][j])
            {
                edge = make_pair(i, j);
                positive.push_back(edge);
            }
        }
    }

    //fout << "\\begin{figure}" << endl
    //<< "\\centering" << endl
    fout << "\\begin{tikzpicture}[scale=.5]" << endl;
    int q = (int)positive.size();
    for(int i = 0; i < q; i++)
    {
        fout << "\\node[vertex] (" << positive[i].first << "_" << positive[i].second
        << ") at (" << 90 + i*360/q << ":3) {" << positive[i].first
        << ", " << positive[i].second << "};" << endl;
    }
    fout << endl;

    fout << "\\node at (270:4) {$\\mathcal{Q}$ for $y = ";
    printVec(p, fout);
    fout << "$};" << endl;

    for(int i = 0; i < q; i++)
    {
        pair<int, int> empty;
        vector<pair<int, int> > adj((int)s[positive[i].first].size() + (int)s[positive[i].second].size(), empty);

        merge(s[positive[i].first].begin(), s[positive[i].first].end(), s[positive[i].second].begin(),
            s[positive[i].second].end(), adj.begin(), pairComp);

        vector<pair<int, int> >::iterator it;
        it = unique(adj.begin(), adj.end());
        adj.resize(it-adj.begin());

        for(int j = 0; j < (int)adj.size(); j++)
        {
            if(pairComp(positive[i], adj[j]))
               {
                    fout << "\\draw[positive] (" << positive[i].first << "_" << positive[i].second
                    << ")--(" << adj[j].first << "_" << adj[j].second << ");" << endl;
               }

        }
    }



    fout << "\\end{tikzpicture}" << endl;
    //<< "\\caption{$\\mathcal{Q}$ graph for ";
    //printVec(p, fout);
    //fout << "}" << endl
    //<< "\\end{figure}" << endl << endl;
}

void printDoubleVec(vector<vector<int> > vec)
{
    for(int i = 0; i < (int)vec.size(); i++)
    {
        for(int j = 0; j < (int) vec[i].size(); j++)
        {
            cout << vec[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void printVertexEdges(vector<vector<pair<int, int> > > s)
{
    for(int i = 0; i < (int)s.size(); i++)
    {
        cout << "Vertex " << i << ": " << endl;
        for(int j = 0; j < (int)s[i].size(); j++)
        {
            cout << "(" << s[i][j].first << ", " << s[i][j].second << ") ";
        }
        cout << endl;
    }
}

bool pairEq(pair<int, int> a, pair<int, int> b)
{
    if(a.first == b.first && a.second == b.second)
    {
        return true;
    }
    return false;
}

bool pairComp(pair<int, int> a, pair<int, int> b)
{
    if(a.first > b.first)
    {
        return false;
    }
    if(a.first < b.first)
    {
        return true;
    }
    if(a.second < b.second)
    {
        return true;
    }
    return false;
}

int coloring(vector<vector<int> > &adj, vector<int> &colors)
{
    int n = (int)adj.size();
    if(n==0)
    {
        return 0;
    }
    int k = n + 1;
    int numColors = n;
    vector<int> bestColoring;
    vector<int> perm;
    for(int i = 0; i < n; i++)
    {
        perm.push_back(i);
    }
    do
    {
        numColors = greedyColoring(adj, colors, perm, k);
        if(numColors < k)
        {
            k = numColors;
            bestColoring.assign(colors.begin(), colors.end());
        }
    }while(next_permutation(perm.begin(), perm.end()));
    colors.assign(bestColoring.begin(), bestColoring.end());

    return k;
}

int greedyColoring(vector<vector<int> > &adj, vector<int> &colors, vector<int> &perm, int best)
{
    int n = (int)adj.size();
    vector<int> allColors;
    for(int i = 0; i < (int)adj.size(); i++)
    {
        allColors.push_back(i);
    }
    vector<int> colorsLeft;
    for(int i = 0; i < n; i++)
    {
        colorsLeft = allColors;
        for(int j = 0; j < i; j++)
        {
            if(adj[perm[i]][perm[j]])
            {
                colorsLeft[colors[perm[j]]] = n; //can't use the color cause adjacent
            }
        }
        for(int j = 0; j < n; j++)
        {
            if(colorsLeft[j] < n)
            {
                colors[perm[i]] = colorsLeft[j];
                break;
            }
        }
        if(*max_element(colors.begin(), colors.end()) + 1 > best)
        return best;
    }
    return *max_element(colors.begin(), colors.end()) + 1;
}

pair<vector<vector<int> > , vector<pair<int, int> > > printQColoring(vector<int> & p, int &red)
{
    int n = (int)p.size();

    vector<int> temp(n, 0);
    vector<vector<int> > aPlus(n, temp);
    vector<vector<int> > aMinus = aPlus;
    int numRed = 0;
    for(int i = 1; i < n; i++)
    {
        for(int j = 1; j <= p[i]; j++)
        {
            if(j + i <= n)
            {
                aPlus[j - 1][j + i - 1] = 1;
                aPlus[j + i - 1][j - 1] = 1;
            }
        }
        if(p[i] + i + 1 <= n)
        {
            aMinus[p[i]][p[i] + i] = 1;
            aMinus[p[i] + i][p[i]] = 1;
            numRed++;
        }
    }
    red = numRed;
    //printDoubleVec(aPlus);
    //printDoubleVec(aMinus);

    vector<pair<int, int> > stuff;
    vector<vector<pair<int, int> > > s(n, stuff);
    pair <int, int> edge;

    for(int i = 0; i < n; i++) //for each vertex in P
    {
        for(int j = 0; j < n; j++)
        {
            if(aMinus[i][j]) //if it's connected by a red edge
            {
                for(int k = 0; k < n; k++)
                {
                    if(aPlus[j][k])
                    {
                        if(j < k)
                        {
                            edge = make_pair(j, k);
                        }
                        else
                        {
                            edge = make_pair(k, j);
                        }
                        s[i].push_back(edge);
                    }
                }
            }
        }
        sort(s[i].begin(), s[i].end(), pairComp);
        vector<pair<int, int> >::iterator it;
        it = unique(s[i].begin(), s[i].end(), pairEq);
        s[i].resize(it - s[i].begin());
    }
    //printVertexEdges(s);

    vector<pair<int, int> > positive;
    for(int i = 0; i < n; i++)
    {
        for(int j = i; j < n; j++)
        {
            if(aPlus[i][j])
            {
                edge = make_pair(i, j);
                positive.push_back(edge);
            }
        }
    }

    int q = (int)positive.size();
    vector<int> zero(q, 0);
    vector<vector<int> > adjacency(q, zero);

    for(int i = 0; i < q; i++)
    {
        pair<int, int> empty;
        vector<pair<int, int> > adj((int)s[positive[i].first].size() + (int)s[positive[i].second].size(), empty);

        merge(s[positive[i].first].begin(), s[positive[i].first].end(), s[positive[i].second].begin(),
            s[positive[i].second].end(), adj.begin(), pairComp);

        vector<pair<int, int> >::iterator it;
        it = unique(adj.begin(), adj.end());
        adj.resize(it-adj.begin());

        for(int j = 0; j < (int)adj.size(); j++)
        {
            if(pairComp(positive[i], adj[j]))
               {
                    adjacency[i][index(adj[j], positive)] = 1;
                    adjacency[index(adj[j], positive)][i] = 1;
                    //TODO convert this to an entry in adj and give it to your coloring function
               }

        }
    }
    pair<vector<vector<int> > , vector<pair<int, int> > > result;
    result.first = adjacency;
    result.second = positive;
    return result;
}

int index(pair<int, int> p, vector<pair<int, int> > &positive)
{
    for(int i = 0; i < (int)positive.size(); i++)
    {
        if(pairEq(p, positive[i]))
        {
            return i;
        }
    }
    return (int)positive.size();
}

vector<vector<int> > indetString(vector<int> p, vector<int> &colors, vector<pair<int, int> > &positive, int &k)
{
    int n = (int)p.size();
    int q = (int)positive.size();
    vector<int> emp;
    vector<vector<int> > w(n, emp);
    for(int i = 0; i < q; i++)
    {
        w[positive[i].first].push_back(colors[i]);
        w[positive[i].second].push_back(colors[i]);
    }
    for(int i = 0; i < n; i++)
    {
        sort(w[i].begin(), w[i].end());
        w[i].resize(unique(w[i].begin(), w[i].end()) - w[i].begin());
    }

    vector<int> temp(n, 0);
    vector<vector<int> > aMinus(n, temp);

    for(int i = 1; i < n; i++)
    {
        if(p[i] + i + 1 <= n)
        {
            aMinus[p[i]][p[i] + i] = 1;
            aMinus[p[i] + i][p[i]] = 1;
        }
    }
    vector<int> colorsLeft;
    vector<int> totColors;
    if(!colors.empty())
    {
        for(int i = 0; i < k; i++)
        {
            totColors.push_back(i);
        }
    }

    for(int i = 0; i < n; i++)
    {
        if(w[i].empty())
        {
            colorsLeft = totColors;
            for(int j = 0; j < n; j++)
            {
                if(aMinus[i][j])
                {
                    vector<int>::iterator it = set_difference(colorsLeft.begin(), colorsLeft.end(),
                                                               w[j].begin(), w[j].end(), colorsLeft.begin());
                    colorsLeft.resize(it - colorsLeft.begin());
                }
            }
            if(colorsLeft.empty())
            {
                w[i].push_back((int)totColors.size());
                totColors.push_back((int)totColors.size());
            }
            else
            {
                w[i].push_back(colorsLeft[0]);
            }
        }
    }
    //printIndetString(w);
    k = (int)totColors.size();
    return w;
}

void printIndetString(vector<vector<int> > &w)
{
    for(int i = 0; i < (int)w.size(); i++)
    {
        cout << i << ":  ";
        printVec(w[i], cout);
        cout << endl;
    }
}

