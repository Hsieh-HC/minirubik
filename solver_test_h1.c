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
    STATES = PERMUTATIONS * ORIENTATIONS,
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

static state_t apply_move(state_t state, uint8_t move)
{
    uint8_t face = (uint8_t)(move / 3U);
    uint8_t turns = (uint8_t)(move % 3U + 1U);
    for (uint8_t t = 0; t < turns; ++t) {
        state_t p_part = quarter_turn_p(state, face);
        state_t o_part = quarter_turn_o(state, face);
        state_t next;
        for (uint8_t i = 0; i < CUBIES; ++i) {
            next.p[i] = p_part.p[i];
            next.o[i] = o_part.o[i];
        }
        state = next;
    }
    return state;
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

static int valid(const state_t *state)
{
    uint8_t sum = 0;
    for (uint8_t i = 0; i < CUBIES; ++i) {
        if (state->p[i] >= CUBIES || state->o[i] >= 3)
            return 0;
        for (uint8_t j = 0; j < i; ++j)
            if (state->p[j] == state->p[i])
                return 0;
        sum = (uint8_t) (sum + state->o[i]);
    }
    return sum % 3U == 0;
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
    free(queue);
    if (tail != STATES_O) {
        free(dist);
        return NULL;
    }
    return dist;
}

static int parse_state(const char *input, state_t *state)
{
    for (int i = 0; i < 14; ++i) {
        int limit = i < 7 ? 7 : 3;
        if (input[i] < '1' || input[i] > '0' + limit)
            return 0;
        (i < 7 ? state->p : state->o)[i % 7] = (uint8_t) (input[i] - '1');
    }
    return input[14] == '\0' && valid(state);
}

static int output_failed(void)
{
    return fflush(stdout) != 0 || ferror(stdout);
}

static int self_test(uint8_t *table_p, uint8_t *table_o)
{
    const state_t solved = {{0, 1, 2, 3, 4, 5, 6}, {0}};
    state_t state;
    uint8_t max_p = 0, max_o = 0;
    for (uint8_t move = 0; move < MOVES; ++move) {
        state = solved;
        state = apply_move(state, move);
        state = apply_move(state, inverse_move[move]);
        if (memcmp(&solved, &state, sizeof solved))
            return 0;
    }
    for (uint32_t rank = 0; rank < PERMUTATIONS; ++rank) {
        if (table_p[rank] == UINT8_MAX || table_p[rank] > 7)
            return 0;
        max_p = table_p[rank] > max_p ? table_p[rank] : max_p;
        unrank_permutation(rank, &state);
        if (!valid(&state) || permutation_rank(&state) != rank)
            return 0;
    }
    for (uint32_t rank = 0; rank < ORIENTATIONS; ++rank) {
        if (table_o[rank] == UINT8_MAX || table_o[rank] > 6)
            return 0;
        unrank_orientation(rank, &state);
        max_o = table_o[rank] > max_o ? table_o[rank] : max_o;
        if (!valid(&state) || orientation_rank(&state) != rank)
            return 0;
    }
    if (table_p[0] != 0 || table_o[0] != 0 || max_p != 7 || max_o != 6)
        return 0;
    return 1;
}

static uint8_t heuristic(uint8_t *table_p, uint8_t *table_o, state_t state){
    uint16_t p_rank = (uint16_t) permutation_rank(&state);
    uint16_t o_rank = (uint16_t) orientation_rank(&state);
    uint8_t hp = table_p[p_rank];
    uint8_t ho = table_o[o_rank];
    uint8_t h = hp >= ho ? hp : ho;
    return h;
}

static uint8_t *build_exact_distance_table(void)
{
    uint8_t *dist = malloc(STATES);
    uint32_t *queue = malloc((size_t)STATES * sizeof *queue);

    if (!dist || !queue) {
        free(dist);
        free(queue);
        return NULL;
    }

    uint16_t permutation[3][PERMUTATIONS];
    uint16_t orientation[3][ORIENTATIONS];

    state_t state;

    /* Build permutation quarter-turn transitions. */
    for (uint16_t rank = 0; rank < PERMUTATIONS; ++rank) {
        unrank_permutation(rank, &state);

        for (uint8_t face = 0; face < 3; ++face) {
            state_t next = quarter_turn_p(state, face);
            permutation[face][rank] =
                (uint16_t)permutation_rank(&next);
        }
    }

    /* Build orientation quarter-turn transitions. */
    for (uint16_t rank = 0; rank < ORIENTATIONS; ++rank) {
        unrank_orientation(rank, &state);

        for (uint8_t face = 0; face < 3; ++face) {
            state_t next = quarter_turn_o(state, face);
            orientation[face][rank] =
                (uint16_t)orientation_rank(&next);
        }
    }

    memset(dist, UINT8_MAX, STATES);

    uint32_t head = 0;
    uint32_t tail = 1;

    dist[0] = 0;
    queue[0] = 0;

    while (head < tail) {
        uint32_t here = queue[head++];

        uint16_t p = (uint16_t)(here / ORIENTATIONS);
        uint16_t o = (uint16_t)(here % ORIENTATIONS);

        for (uint8_t face = 0; face < 3; ++face) {
            uint16_t next_p = p;
            uint16_t next_o = o;

            for (uint8_t turn = 0; turn < 3; ++turn) {
                next_p = permutation[face][next_p];
                next_o = orientation[face][next_o];

                uint32_t there =
                    (uint32_t)next_p * ORIENTATIONS + next_o;

                if (dist[there] == UINT8_MAX) {
                    dist[there] = (uint8_t)(dist[here] + 1);
                    queue[tail++] = there;
                }
            }
        }
    }

    free(queue);

    if (tail != STATES) {
        free(dist);
        return NULL;
    }

    return dist;
}

static int test_h1(const uint8_t *table_p,
                   const uint8_t *table_o,
                   const uint8_t *exact_dist)
{
    uint32_t checked = 0;

    for (uint32_t p = 0; p < PERMUTATIONS; ++p) {
        for (uint32_t o = 0; o < ORIENTATIONS; ++o) {

            uint32_t full_rank = p * ORIENTATIONS + o;

            uint8_t hp = table_p[p];
            uint8_t ho = table_o[o];
            uint8_t h = hp >= ho ? hp : ho;

            uint8_t d = exact_dist[full_rank];

            if (h > d) {
                fprintf(stderr,
                        "H1 failed: rank=%u p=%u o=%u h=%u d=%u\n",
                        full_rank, p, o, h, d);
                return 0;
            }

            ++checked;
        }
    }

    printf("H1 passed: %u states checked\n", checked);
    return 1;
}

int main(int argc, char **argv)
{
    state_t start;
    uint8_t diameter_p, diameter_o;

    if (argc == 2 && !strcmp(argv[1], "--h1-test")) {
        uint8_t diameter_p, diameter_o;

        uint8_t *table_p = build_table_p(&diameter_p);
        uint8_t *table_o = build_table_o(&diameter_o);
        uint8_t *exact_dist = build_exact_distance_table();

        if (!table_p || !table_o || !exact_dist) {
            free(table_p);
            free(table_o);
            free(exact_dist);
            fputs("could not build H1 tables\n", stderr);
            return 1;
        }

        int ok = test_h1(table_p, table_o, exact_dist);

        free(table_p);
        free(table_o);
        free(exact_dist);

        return ok ? 0 : 1;
    }

    if (argc == 2 && !strcmp(argv[1], "--self-test")) {
        uint8_t *table_p = build_table_p(&diameter_p);
        uint8_t *table_o = build_table_o(&diameter_o);
        if (!table_p || !table_o) {
            fputs("could not build heuristic tables\n", stderr);
            return 1;
        }
        if (!self_test(table_p, table_o)) {
            fputs("self-test failed\n", stderr);
            return 1;
        }
        free(table_p);
        free(table_o);
        if (diameter_p != 7 || diameter_o != 6) {
            fputs("BFS check failed\n", stderr);
            return 1;
        }
        puts("5040 permutations; diameter 7");
        puts("729 orientations; diameter 6");
        return output_failed();
    }
    if (argc != 2 || !parse_state(argv[1], &start)) {
        fprintf(stderr, "usage: %s PPPPPPPOOOOOOO\n",
                argc > 0 && argv[0] ? argv[0] : "solver");
        return 2;
    }

    /* 2. build heuristic tables */
    uint8_t *table_p = build_table_p(&diameter_p);
    uint8_t *table_o = build_table_o(&diameter_o);
    if (!table_p || !table_o) {
        free(table_p);
        free(table_o);
        fputs("could not build heuristic tables\n", stderr);
        return 1;
    }

    /* 3. first IDA* threshold */
    uint8_t threshold = heuristic(table_p, table_o, start);

    /* iterative IDA* */
    uint8_t path[12];
    state_t states[12];
    uint8_t next_move[12];
    int depth;
    uint8_t next_threshold;
    int found = 0;
    while (!found) {
        depth = 0;
        states[0] = start;
        next_move[0] = 0;
        next_threshold = 20;
        while (depth >= 0) {
            state_t current = states[depth];
            uint16_t p_rank = (uint16_t) permutation_rank(&current);
            uint16_t o_rank = (uint16_t) orientation_rank(&current);
            uint8_t hp = table_p[p_rank];
            uint8_t ho = table_o[o_rank];
            uint8_t h = hp >= ho ? hp : ho;
            uint8_t f = depth + h;
            if (f > threshold) {
                next_threshold = next_threshold < f ? next_threshold : f;
                depth--;
                continue;
            }
            if (h == 0) {
                found = 1;
                break;
            }
            if (next_move[depth] >= MOVES) {
                depth--;
                continue;
            }
            uint8_t move = next_move[depth]++;
            uint8_t face = move / 3U;
            if (depth > 0) {
                uint8_t last_face = path[depth - 1] / 3U;
                if (face == last_face) {
                    continue; 
                }
            }
            states[depth + 1] = apply_move(current, move);
            path[depth] = move;
            next_move[depth + 1] = 0;
            depth++;
        }
        if (!found) {
            threshold = next_threshold;
        }
    }
    for (int i = 0; i < depth; i++) {
        if (i > 0) {
            putchar(' ');
        }
        printf("%s", move_names[path[i]]);
    }
    putchar('\n');
    state_t check = start;
    for (int i = 0; i < depth; i++) {
        check = apply_move(check, path[i]);
    }
    if (permutation_rank(&check) != 0 || orientation_rank(&check) != 0) {
        fputs("solution verification failed\n", stderr);
        free(table_p);
        free(table_o);
        return 1;
    }
    free(table_p);
    free(table_o);
    return output_failed();
}