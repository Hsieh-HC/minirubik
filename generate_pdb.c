#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum {
    CUBIES = 7,
    PERMUTATIONS = 5040,
    STATES_P = PERMUTATIONS,
    ORIENTATIONS = 729,
    STATES_O = ORIENTATIONS,
    MOVES = 9
};

typedef struct {
    uint8_t p[CUBIES], o[CUBIES];
} state_t;

static const char *const move_names[MOVES] = {"R",  "R2", "R'", "B", "B2",
                                              "B'", "D",  "D2", "D'"};
static const uint8_t inverse_move[MOVES] = {2, 1, 0, 5, 4, 3, 8, 7, 6};
static const uint8_t source[3][CUBIES] = {
    {1, 4, 2, 0, 3, 5, 6},
    {0, 1, 2, 4, 5, 6, 3},
    {0, 2, 5, 3, 1, 4, 6},
};

static const uint8_t twist[3][CUBIES] = {
    {1, 2, 0, 2, 1, 0, 0},
    {0, 0, 0, 1, 2, 1, 2},
    {0, 0, 0, 0, 0, 0, 0},
};

static state_t quarter_turn_p(state_t state, uint8_t face)
{
    state_t result;

    for (uint8_t i = 0; i < CUBIES; ++i) {
        uint8_t from = source[face][i];
        result.p[i] = state.p[from];
    }
    return result;
}
static state_t quarter_turn_o(state_t state, uint8_t face)
{
    state_t result;

    for (uint8_t i = 0; i < CUBIES; ++i) {
        uint8_t from = source[face][i];
        result.o[i] = (uint8_t) ((state.o[from] + twist[face][i]) % 3U);
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
static uint32_t orientation_rank(const state_t *state)
{
    uint32_t o = 0;
    for (uint8_t i = 0; i < 6; ++i)
        o = o * 3U + state->o[i];
    return o;
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
        state->o[i] = 0;
    }
}
static void unrank_orientation(uint32_t rank, state_t *state)
{
    uint32_t o = rank;
    uint8_t sum = 0;
    for (uint8_t i = 6; i-- > 0;) {
        state->o[i] = (uint8_t) (o % 3U);
        sum = (uint8_t) (sum + state->o[i]);
        o /= 3U;
        state->p[i] = (uint8_t) i;
    }
    state->o[6] = (uint8_t) ((3U - sum % 3U) % 3U);
    state->p[6] = 6;
}

static uint8_t *build_table_p(uint8_t *diameter)
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
            state_t next = quarter_turn_p(state, face);
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
    printf("static const uint8_t table_p[%d] = {\n", PERMUTATIONS);

    for (int i = 0; i < PERMUTATIONS; ++i) {
        printf("%u", dist[i]);

        if (i + 1 != PERMUTATIONS)
            printf(", ");

        if ((i + 1) % 16 == 0)
            printf("\n");
    }

    printf("\n};\n");
    free(queue);
    if (tail != STATES_P) {
        free(dist);
        return NULL;
    }
    return dist;
}
static uint8_t *build_table_o(uint8_t *diameter)
{
    uint8_t *dist = malloc(STATES_O);
    uint32_t *queue = malloc((size_t) STATES_O * sizeof *queue);
    if (!dist || !queue) {
        free(dist);
        free(queue);
        return NULL;
    }
    uint16_t orientation[3][ORIENTATIONS];
    uint32_t head = 0, tail = 1, level_end = 1;
    state_t state;
    for (uint16_t rank = 0; rank < ORIENTATIONS; ++rank) {
        unrank_orientation(rank, &state);
        for (uint8_t face = 0; face < 3; ++face) {
            state_t next = quarter_turn_o(state, face);
            orientation[face][rank] =
                (uint16_t) (orientation_rank(&next) % ORIENTATIONS);
        }
    }
    memset(dist, UINT8_MAX, STATES_O);
    dist[0] = 0;
    queue[0] = 0;
    *diameter = 0;
    while (head < tail){
        if (head == level_end) {
            level_end = tail;
            ++*diameter;
        }
        uint32_t here = queue[head++];
        uint16_t o = (uint16_t) here;
        for (uint8_t face = 0; face < 3; ++face) {
            uint16_t next_o = o;
            for (uint8_t turn = 0; turn < 3; ++turn) {
                next_o = orientation[face][next_o];
                uint32_t there = (uint32_t) next_o;
                if (dist[there] == UINT8_MAX) {
                    dist[there] = dist[here] + 1;
                    queue[tail++] = there;
                }
            }
        }
    }
    printf("static const uint8_t table_o[%d] = {\n", ORIENTATIONS);

    for (int i = 0; i < ORIENTATIONS; ++i) {
        printf("%u", dist[i]);

        if (i + 1 != ORIENTATIONS)
            printf(", ");

        if ((i + 1) % 16 == 0)
            printf("\n");
    }

    printf("\n};\n");
    free(queue);
    if (tail != STATES_O) {
        free(dist);
        return NULL;
    }
    return dist;
}


int main(void)
{
    uint8_t diameter_p, diameter_o;

    uint8_t *p = build_table_p(&diameter_p);
    uint8_t *o = build_table_o(&diameter_o);

    if (!p || !o) {
        free(p);
        free(o);
        return 1;
    }

    free(p);
    free(o);
    return 0;
}