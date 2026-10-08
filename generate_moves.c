#include <stdio.h>
#include <stdint.h>

#define FACES 3
#define CORNERS 7

static const uint8_t source[FACES][CORNERS] = {
    {1, 4, 2, 0, 3, 5, 6},
    {0, 1, 2, 4, 5, 6, 3},
    {0, 2, 5, 3, 1, 4, 6},
};

static const uint8_t twist[FACES][CORNERS] = {
    {1, 2, 0, 2, 1, 0, 0},
    {0, 0, 0, 1, 2, 1, 2},
    {0, 0, 0, 0, 0, 0, 0},
};

int main(void)
{
    for (int face = 0; face < 3; face++) {
        uint8_t cur_source[7] = {0,1,2,3,4,5,6};
        uint8_t cur_twist[7]  = {0,0,0,0,0,0,0};

        for (int turns = 1; turns <= 3; turns++) {
            uint8_t next_source[7];
            uint8_t next_twist[7];

            for (int i = 0; i < 7; i++) {
                int from = source[face][i];

                next_source[i] = cur_source[from];
                next_twist[i] =
                    (cur_twist[from] + twist[face][i]) % 3;
            }

            for (int i = 0; i < 7; i++) {
                cur_source[i] = next_source[i];
                cur_twist[i] = next_twist[i];
            }

            printf("move %d source: ", face * 3 + turns - 1);
            for (int i = 0; i < 7; i++)
                printf("%d%s", cur_source[i],
                       i == 6 ? "\n" : ",");

            printf("move %d twist : ", face * 3 + turns - 1);
            for (int i = 0; i < 7; i++)
                printf("%d%s", cur_twist[i],
                       i == 6 ? "\n" : ",");
        }
    }

    return 0;
}