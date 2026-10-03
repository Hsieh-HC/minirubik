#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum {
    CUBIES = 7,
    ORIENTATIONS = 729,
    STATES_O = ORIENTATIONS,
};

typedef struct {
    uint8_t o[CUBIES];
} state_t;


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
static state_t quarter_turn(state_t state, uint8_t face)
{
    state_t result;

    for (uint8_t i = 0; i < CUBIES; ++i) {
        uint8_t from = source[face][i];
        result.o[i] = (uint8_t) ((state.o[from] + twist[face][i]) % 3U);
    }
    return result;
}
static uint32_t orientation_rank(const state_t *state)
{
    uint32_t o = 0;
    for (uint8_t i = 0; i < 6; ++i)
        o = o * 3U + state->o[i];
    return o;
}
static void unrank_orientation(uint32_t rank, state_t *state)
{
    uint32_t o = rank;
    uint8_t sum = 0;
    for (uint8_t i = 6; i-- > 0;) {
        state->o[i] = (uint8_t) (o % 3U);
        sum = (uint8_t) (sum + state->o[i]);
        o /= 3U;
    }
    state->o[6] = (uint8_t) ((3U - sum % 3U) % 3U);
}
static uint8_t *build_table(uint8_t *diameter)
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
            state_t next = quarter_turn(state, face);
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
    free(queue);
    if (tail != STATES_O) {
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
        fputs("could not build orientation table\n", stderr);
        return 1;
    }
    printf("%d, %d\n", STATES_O, diameter);
    free(table);
    return 0;
}