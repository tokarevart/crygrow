#include <iostream>
#include <unordered_map>
#include <memory>
#include <functional>

struct base {
    int a = 0;

    virtual void foo() = 0;
};

struct derived0 : base {
    int b = 1;

    void foo() override {
        std::cout << "derived0::foo" << std::endl;
    }

    void bar() {
        std::cout << "derived0::bar" << std::endl;
    }

    void foo_call() {
        derived0::foo();
    }
};

struct derived1 : derived0 {
    int b = 1;
    int c = 2;
    
    void foo() override {
        std::cout << "derived1::foo" << std::endl;
    }

    void bar() {
        std::cout << "derived1::bar" << std::endl;
    }

    derived1() {
        derived0::b++;
        std::cout << "derived1::b " << b << std::endl;
        std::cout << "derived0::b " << derived0::b << std::endl;

        derived0::a++;
        std::cout << "derived0::a " << derived0::a << std::endl;
        std::cout << "base::a " << base::a << std::endl;
    }
};

class nexter {
public:
    std::function<int(int)> get() {
        return [this](int val) { return next(val); };
    }
    
    nexter(int shift) : m_shift{shift} {}


private:
    int m_shift;

    int next(int val) {
        return val + m_shift;
    }
};

template <typename T>
struct cls {
    int a = 1;
    float b = 0.5;

    auto f() const {
        if constexpr (std::is_same_v<kek, T>)
            return a;
        else
            return b;
    }

    cls() {}
};


int main() {
    //derived1* ptr = new derived1;
    //ptr->derived0::a = 2;
    //ptr->foo();
    //ptr->foo_call();
    //std::cout << dynamic_cast<derived1*>(ptr)->derived0::a << ' ' << dynamic_cast<derived1*>(ptr)->b;

    //std::unordered_map<int, std::unique_ptr<int>> example;
    //example.emplace(1, std::make_unique<int>(1));
    //example.emplace(2, std::make_unique<int>(2));
    //auto search = example.find(2);
    //if (search != example.end())
    //    std::cout << "Found " << search->first << " " << search->second << '\n';
    //else
    //    std::cout << "Not found\n";

    //auto next = nexter(2).get();
    //std::cout << next(2) << ' ' << next(next(3));

    std::cout << cls<int>().f() << ' ' << cls<float>().f();

    return 0;
}
