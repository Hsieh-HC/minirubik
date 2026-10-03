#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum {
    CUBIES = 7,
    PERMUTATIONS = 5040,
    STATES_P = PERMUTATIONS,
};

typedef struct {
    uint8_t p[CUBIES];
} state_t;


static const uint8_t source[3][CUBIES] = {
    {1, 4, 2, 0, 3, 5, 6},
    {0, 1, 2, 4, 5, 6, 3},
    {0, 2, 5, 3, 1, 4, 6},
};
static state_t quarter_turn(state_t state, uint8_t face)
{
    state_t result;

    for (uint8_t i = 0; i < CUBIES; ++i) {
        uint8_t from = source[face][i];
        result.p[i] = state.p[from];
    }
    return result;
}
static uint32_t permutation_rank(const state_t *state)
{
    uint32_t p = 0;
    for (uint8_t i = 0; i < CUBIES; ++i) {
        uint8_t smaller = 0;
        for (uint8_t j = (uint8_t) (i + 1U); j < CUBIES; ++j)
            if (state->p[j] < state->p[i])
                ++smaller;
        p = p * (CUBIES - i) + smaller;
    }
    return p;
}
static void unrank_permutation(uint32_t rank, state_t *state)
{
    uint8_t available[CUBIES] = {0, 1, 2, 3, 4, 5, 6};
    uint32_t p = rank, f = 720;
    for (uint8_t i = 0; i < CUBIES; ++i) {
        uint8_t q = (uint8_t) (p / f);
        p %= f;
        state->p[i] = available[q];
        for (uint8_t j = q; j + 1U < CUBIES - i; ++j)
            available[j] = available[j + 1U];
        if (i < 5)
            f /= 6U - i;
    }
}
static uint8_t *build_table(uint8_t *diameter)
{
    uint8_t *dist = malloc(STATES_P);
    uint32_t *queue = malloc((size_t) STATES_P * sizeof *queue);
    if (!dist || !queue) {
        free(dist);
        free(queue);
        return NULL;
    }
    uint16_t permutation[3][PERMUTATIONS];
    uint32_t head = 0, tail = 1, level_end = 1;
    state_t state;
    for (uint16_t rank = 0; rank < PERMUTATIONS; ++rank) {
        unrank_permutation((uint32_t) rank, &state);
        for (uint8_t face = 0; face < 3; ++face) {
            state_t next = quarter_turn(state, face);
            permutation[face][rank] =
                (uint16_t) (permutation_rank(&next));
        }
    }
    memset(dist, UINT8_MAX, STATES_P);
    dist[0] = 0;
    queue[0] = 0;
    *diameter = 0;
    while (head < tail){
        if (head == level_end) {
            level_end = tail;
            ++*diameter;
        }
        uint32_t here = queue[head++];
        uint16_t p = (uint16_t) here;
        for (uint8_t face = 0; face < 3; ++face) {
            uint16_t next_p = p;
            for (uint8_t turn = 0; turn < 3; ++turn) {
                next_p = permutation[face][next_p];
                uint32_t there = (uint32_t) next_p;
                if (dist[there] == UINT8_MAX) {
                    dist[there] = dist[here] + 1;
                    queue[tail++] = there;
                }
            }
        }
    }
    free(queue);
    if (tail != STATES_P) {
        free(dist);
        return NULL;
    }
    return dist;
}
int main(void)
{
    uint8_t diameter;
    uint8_t *table = build_table(&diameter);
    if (!table){
        fputs("could not build permutation table\n", stderr);
        return 1;
    }
    printf("%d, %d\n", STATES_P, diameter);
    free(table);
    return 0;
}