
    int main() {
        #if defined(_M_ARM)
            return 2;
        #elif defined(_M_ARM64)
            return 3;
        #elif defined(_M_AMD64)
            return 4;
        #elif defined(_M_X64)
            return 5;
        #elif defined(_M_IX86)
            return 6;
        #else
            return 0;
    #endif
}
