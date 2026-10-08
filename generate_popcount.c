#include <stdio.h>

int main(void)
{
    printf("popcount_table:\n");

    for (int i = 0; i < 128; i++) {
        int x = i;
        int count = 0;

        while (x != 0) {
            count += x & 1;
            x >>= 1;
        }

        if (i % 16 == 0)
            printf("    .byte ");

        printf("%d", count);

        if (i % 16 == 15)
            printf("\n");
        else
            printf(",");
    }

    return 0;
}