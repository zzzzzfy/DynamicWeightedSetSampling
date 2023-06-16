#include<cstdio>
#include<iostream>
#include<vector>
#include<bitset>
#include "size_block.hpp"
#include "basic_alias.hpp"
#ifndef DYNAMICSETSAMPLING_CHUNKSCAPEGOATTREE_HPP
#define DYNAMICSETSAMPLING_CHUNKSCAPEGOATTREE_HPP
bool use_alias;

template<typename Method>
class ChunkScapegoatTree {

    struct chunk {
        aliasMethod* alias;
        std::vector<Element> elements;
        int size;
        float totWeight;
        float U;
        int L, R;

        chunk(std::vector<Element>& _elements) {
            elements = _elements;
            size = elements.size();
            L = R = elements[0].key;
            totWeight = 0;
            for (auto e : elements)
            {
                totWeight += e.weight;
                L = min(L, e.key);
                R = max(R, e.key);
            }
            if(use_alias)
            alias = new aliasMethod(size, elements);
            sort(elements.begin(), elements.end(), cmp_element_key);
        }

        chunk(std::vector<Element>& _elements, float _U) {
            elements = _elements;
            size = elements.size();
            U = _U;
            L = R = elements[0].key;
            totWeight = 0;
            for (auto e : elements)
            {
                totWeight += e.weight;
                L = min(L, e.key);
                R = max(R, e.key);
            }
            if (use_alias)
            alias = new aliasMethod(size, elements);
            sort(elements.begin(), elements.end(), cmp_element_key);
        }

        ~chunk() {
            elements.clear();
            delete(alias);
        }
        void insert(Element element) {
            L = min(element.key, L);
            R = max(element.key, R);
            totWeight += element.weight;
            elements.insert(upper_bound(elements.begin(), elements.end(), element, cmp_element_key_struct()), element);
            if (use_alias)
            alias->insert_element(element);
            size++;
        }
        void erase(int key) {
            for (auto it = elements.begin(); it != elements.end(); it++)
                if ((*it).key == key) {
                    size--;
                    totWeight -= (*it).weight;
                    elements.erase(it);
                    break;
                }
            if (use_alias)
            alias->delete_element(key);
        }
        int randomSampleValue(mt19937& rng) {
            //
            //return alias->random_sample_value();
            uniform_real_distribution<float> db(0, totWeight);
            float sampl = db(rng);
            for (int i = 0; i < size; i++)
                if (elements[i].weight >= sampl)return elements[i].value;
                else sampl -= elements[i].weight;
        }

        void getElements(int l, int r, vector<Element>& retElements) {
            for (const Element& element : elements)
                if (l <= element.key && element.key <= r)
                    retElements.push_back(element);
        }
    };

    struct node {
        Method* method;
        int L, R, size;
        vector<Element> elements;
        node* leftSon, * rightSon;
        ~node() {
            elements.clear();
            delete(method);
        }
    };
public:
    void setChunkSize(int _cSize) {
        chunkSize = _cSize;
        if (_cSize < 1000) use_alias = false;
        else use_alias = true;
    }
    void init(Element* dataSet, int num) {
        //if(num > 2e9)
        root = buildTree(dataSet, num);
    }
    void init(Element* dataSet, int num, float u) {
        U = u;
        root = buildTree(dataSet, num);
    }
    void insert(Element newData) {
        sta.clear();
        node* nowNode = root, * father = nullptr;
        while (nowNode->leftSon != nullptr)
        {
            father = nowNode;
            sta.push_back(nowNode);
            nowNode->L = min(nowNode->L, newData.key);
            nowNode->R = max(nowNode->R, newData.key);
            if (newData.key < nowNode->rightSon->L)
                nowNode = nowNode->leftSon;
            else nowNode = nowNode->rightSon;
        }
        deleteChunk(nowNode->L);
        nowNode->L = min(nowNode->L, newData.key);
        nowNode->R = max(nowNode->R, newData.key);
        int chunkId = nowNode->elements[0].value;
        chunk* nowChunk = chunks[chunkId];
        nowChunk->insert(newData);
        if (nowChunk->size >= 2 * chunkSize) {
            split(father, nowNode, nowChunk);
            balance();
        }
        else {
            nowNode->elements[0] = chunkToElement(nowChunk, chunkId);
            insertChunk(nowNode->elements[0]);
        }
    }
    void erase(int position) {
        node* nowNode = root;
        while (nowNode->leftSon != nullptr)
        {
            if (position < nowNode->rightSon->L)
                nowNode = nowNode->leftSon;
            else nowNode = nowNode->rightSon;
        }
        deleteChunk(nowNode->L);
        int chunkId = nowNode->elements[0].value;
        chunk* nowChunk = chunks[chunkId];
        nowChunk->erase(position);
        nowNode->elements[0] = chunkToElement(nowChunk, chunkId);
        insertChunk(nowNode->elements[0]);
    }
    vector<int> ask(int l, int r, int k) {
        vector<Method*> sampleMethod;
        getInterval(root, l, r, sampleMethod);
        vector<Element> sampleElements;
        for (int id = 0; id < sampleMethod.size(); id++) {
            auto method = sampleMethod[id];
            sampleElements.push_back(Element(-1, id, method->tot_weight));
        }
        int methodNum = sampleElements.size();

        chunk* chunk1 = getChunk(root, l);
        chunk* chunk2 = getChunk(root, r);
        if (chunk1 == chunk2 && chunk1 != nullptr) chunk1->getElements(l, r, sampleElements);
        else {
            if (chunk1 != nullptr) chunk1->getElements(l, r, sampleElements);
            if (chunk2 != nullptr) chunk2->getElements(l, r, sampleElements);
        }
        for (int i = methodNum; i < sampleElements.size(); i++)
            sampleElements[i].value += methodNum;
        aliasMethod alias(sampleElements.size(), sampleElements);

        vector<int>ans;
        vector<int>sampleNode;
        sampleNode.resize(methodNum, 0);
        for (int i = 0; i < k; i++)
        {
            int sampleRes = alias.random_sample_value();
            if (sampleRes >= methodNum) ans.push_back(sampleRes - methodNum);
            else {
                sampleNode[sampleRes]++;
                //                int chunkId = sampleMethod[sampleRes]->random_sample_value();
                //                ans.push_back(chunks[chunkId]->randomSampleValue(rng));
            }
        }
        for (int i = 0; i < methodNum; i++)
            for (int j = 0; j < sampleNode[i]; j++)
            {
                int chunkId = sampleMethod[i]->random_sample_value();
                ans.push_back(chunks[chunkId]->randomSampleValue(rng));
            }

        return ans;
    }

    void printTree()
    {
        printTreeInfo(root);
    }
private:

    node* root = nullptr;
    std::vector<node*> buildBuf;
    std::vector<node*> sta;
    std::vector<chunk*> chunks;

    mt19937 rng;
    const double alpha = 0.7;
    int chunkSize;
    float U = -1;

    Element chunkToElement(chunk* ch, int id) {
        return Element(ch->L, id, ch->totWeight);
    }

    void deleteChunk(int key)
    {
        node* nowNode = root;
        while (nowNode->leftSon != nullptr)
        {
            nowNode->method->delete_element(key);
            nowNode->size--;
            if (key < nowNode->rightSon->L)
                nowNode = nowNode->leftSon;
            else nowNode = nowNode->rightSon;
        }
        nowNode->method->delete_element(key);
    }
    void insertChunk(Element element)
    {
        node* nowNode = root;
        while (nowNode->leftSon != nullptr)
        {
            nowNode->method->insert_element(element);
            nowNode->size++;
            if (element.key < nowNode->rightSon->L)
                nowNode = nowNode->leftSon;
            else nowNode = nowNode->rightSon;
        }
        nowNode->method->insert_element(element);
    }
    void split(node* fa, node* nowNode, chunk* nowChunk) {
        vector<Element> preElements = nowChunk->elements;
        vector<Element> LElements, RElements;
        int num = preElements.size();
        for (int i = 0; i < num; i++)
            if (i + 1 <= num / 2) LElements.push_back(preElements[i]);
            else RElements.push_back(preElements[i]);
        delete(nowChunk);
        chunk* LChunk, * RChunk;
        if (U == -1) {
            LChunk = new chunk(LElements);
            RChunk = new chunk(RElements);
        }
        else {
            LChunk = new chunk(LElements, U);
            RChunk = new chunk(RElements, U);
        }
        int LId = chunks.size();
        int RId = LId + 1;
        chunks.push_back(LChunk);
        chunks.push_back(RChunk);
        insertChunk(chunkToElement(LChunk, LId));
        insertChunk(chunkToElement(RChunk, RId));
        node* ls = newLeaf(LChunk, LId);
        node* rs = newLeaf(RChunk, RId);
        node* nNode = newNode(ls, rs);
        if (fa->leftSon == nowNode)
            fa->leftSon = nNode;
        else fa->rightSon = nNode;
        delete(nowNode);
    }

    void merge(vector<Element>& FElements, vector<Element> const& LElements, vector<Element> RElements) {
        int LPoint = 0, RPoint = 0;
        int LSize = LElements.size(), RSize = RElements.size();
        while (LPoint < LSize && RPoint < RSize)
        {
            if (LElements[LPoint].weight > RElements[RPoint].weight)
                FElements.push_back(LElements[LPoint++]);
            else FElements.push_back((RElements[RPoint++]));
        }
        while (LPoint < LSize) FElements.push_back(LElements[LPoint++]);
        while (RPoint < RSize) FElements.push_back(RElements[RPoint++]);
    }

    void clear(node* node) {
        if (node->leftSon == nullptr) {
            chunk* nowChunk = chunks[node->elements[0].value];
            if (nowChunk->size != 0)
                buildBuf.push_back(node);
            else {
                delete(node);
                delete(nowChunk);
            }
            return;
        }
        clear(node->leftSon);
        clear(node->rightSon);
        delete(node);
    }


    void update(node* nowNode) {
        node* ls = nowNode->leftSon;
        node* rs = nowNode->rightSon;
        nowNode->size = ls->size + rs->size;
        vector<Element> tmp;
        if (ls->size == 0)
            merge(nowNode->elements, tmp, rs->elements);
        else if (rs->size == 0) merge(nowNode->elements, ls->elements, tmp);
        else merge(nowNode->elements, ls->elements, rs->elements);
        if (nowNode->method->method_name() == "size_block")
            nowNode->method = new Method(nowNode->size, nowNode->elements);
        else nowNode->method = new Method(nowNode->size, U, nowNode->elements);
        nowNode->L = ls->L;
        nowNode->R = rs->R;
    }

    node* newNode(node* ls, node* rs) {
        node* nowNode = new node();
        nowNode->leftSon = ls;
        nowNode->rightSon = rs;
        update(nowNode);
        return nowNode;
    }
    node* newLeaf(chunk* nowChunk, int id) {
        Element nowData(nowChunk->L, id, nowChunk->totWeight);
        node* nowNode = new node();
        nowNode->L = nowChunk->L;
        nowNode->R = nowChunk->R;
        nowNode->elements.push_back(nowData);
        if (nowNode->method->method_name() == "size_block")
            nowNode->method = new Method(1, nowNode->elements);
        else nowNode->method = new Method(1, U, nowNode->elements);
        nowNode->size = 1;
        return nowNode;
    }
    node* build(int l, int r) {
        if (l == r) return buildBuf[l];
        int mid = (l + r) / 2;
        node* ls = build(l, mid);
        node* rs = build(mid + 1, r);
        node* nd = newNode(ls, rs);
        if (ls->leftSon != nullptr) ls->elements.clear();
        if (rs->rightSon != nullptr) rs->elements.clear();
        return nd;
    }
    node* buildTree(Element* dataSet, int num) {
        buildBuf.clear();
        int size = 0;
        int index = 0;
        for (int id = 0; id < num; id += chunkSize) {
            size++;
            vector<Element> elements;
            for (int j = id; j < min(id + chunkSize, num); j++)
                elements.push_back(dataSet[j]);
            chunk* nowChunk;
            if (U == -1) nowChunk = new chunk(elements);
            else nowChunk = new chunk(elements, U);
            node* nd = newLeaf(nowChunk, index++);
            buildBuf.push_back(nd);
            chunks.push_back(nowChunk);
        }
        return build(0, size - 1);
        buildBuf.clear();
    }
    node* rebuildSubtree(node* node) {
        buildBuf.clear();
        clear(node);
        int num = buildBuf.size();
        return build(0, num - 1);
    }
    void balance() {
        node* pre = nullptr;
        int stSize = sta.size();
        for (int i = 0; i < stSize; i++)
        {
            node* pNode = sta[i];
            double size = pNode->size;
            int lsize = pNode->leftSon->size;
            int rsize = pNode->rightSon->size;
            if (size * alpha < max(lsize, rsize))
            {
                if (pNode == root)
                    root = rebuildSubtree(pNode);
                else {
                    if (pre->leftSon == pNode)
                        pre->leftSon = rebuildSubtree(pNode);
                    else pre->rightSon = rebuildSubtree(pNode);
                }
                break;
            }
            pre = pNode;
        }
        sta.clear();
    }

    void getInterval(node* nowNode, int l, int r, vector<Method*>& ans) {
        if (l < nowNode->L && nowNode->R < r) {
            ans.push_back(nowNode->method);
            return;
        }
        if (nowNode->leftSon == nullptr) return;
        if (l<nowNode->leftSon->R && r>nowNode->leftSon->L) getInterval(nowNode->leftSon, l, r, ans);
        if (nowNode->rightSon->L < r && l < nowNode->rightSon->R) getInterval(nowNode->rightSon, l, r, ans);
    }

    chunk* getChunk(node* nowNode, int position) {
        if (nowNode->leftSon == nullptr)
        {
            return chunks[nowNode->elements[0].value];
        }
        if (nowNode->leftSon->L < position && position < nowNode->leftSon->R) return getChunk(nowNode->leftSon, position);
        if (nowNode->rightSon->L < position && position < nowNode->rightSon->R) return getChunk(nowNode->rightSon, position);
        return nullptr;
    }

    void printTreeInfo(node* nowNode) {
        if (nowNode->leftSon == nullptr)
        {
            printf("%d %d %.2f\n", nowNode->elements[0].key, nowNode->elements[0].value, nowNode->elements[0].weight);
            return;
        }
        printTreeInfo(nowNode->leftSon);
        printTreeInfo(nowNode->rightSon);
    }


};


#endif //DYNAMICSETSAMPLING_CHUNKSCAPEGOATTREE_HPP