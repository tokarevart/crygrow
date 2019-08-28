#include <iostream>

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

struct derived1 : derived0, base {
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

int main() {
    derived1* ptr = new derived1;
    ptr->derived0::a = 2;
    ptr->foo();
    ptr->foo_call();
    std::cout << dynamic_cast<derived1*>(ptr)->derived0::a << ' ' << dynamic_cast<derived1*>(ptr)->b;
    return 0;
}
