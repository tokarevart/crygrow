#include <iostream>

struct base {
    int a = 0;

    virtual void foo() {
        std::cout << "base::foo" << std::endl;
    }
};

struct derived0 : base {
    int b = 1;

    void foo() override final {
        std::cout << "derived0::foo" << std::endl;
    }

    void bar() {
        std::cout << "derived0::bar" << std::endl;
    }
};

struct derived1 : base, derived0 {
    int b = 1;
    int c = 2;
    
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

int main() {
    base* ptr = new derived1;
    ptr->a = 2;
    std::cout << dynamic_cast<derived1*>(ptr)->derived0::a << ' ' << dynamic_cast<derived1*>(ptr)->b;
    return 0;
}
