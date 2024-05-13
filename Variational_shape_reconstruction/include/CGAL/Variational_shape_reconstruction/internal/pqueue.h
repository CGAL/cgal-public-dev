#ifndef PQUEUE_H
#define PQUEUE_H

#include <queue>

namespace qem {

template<typename T, typename TM>
class Custom_priority_queue : public std::priority_queue<T, std::vector<T>, TM>
{
public:

    void remove(const T& value)
    {
        auto it = std::find(this->c.begin(), this->c.end(), value);
        if(it != this->c.end())
        {
            this->c.erase(it);
            std::make_heap(this->c.begin(), this->c.end(), this->comp);
        }
        
        return;
    }

    void remove(const std::vector<T>& values)
    {
        int count = 0;

        for(int i = 0; i < values.size(); i++)
        {
            auto it = std::find(this->c.begin(), this->c.end(), values[i]);
            if(it != this->c.end()){
                this->c.erase(it);
                count++;
            }
        }

        // std::cout << "Queue remove count: " << count << std::endl;
        if(count > 0)
            std::make_heap(this->c.begin(), this->c.end(), this->comp);

        return;
    }
};




} // namespace qem

#endif
