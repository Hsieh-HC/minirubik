#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

enum {
    CUBIES = 7,
    PERMUTATIONS = 5040,
    ORIENTATIONS = 729,
    MOVES = 9
};

typedef struct {
    uint8_t p[CUBIES];
    uint8_t o[CUBIES];
} state_t;

/*
 * Quarter-turn definitions.
 *
 * face 0 = R
 * face 1 = B
 * face 2 = D
 */
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

static const uint8_t face_table[MOVES] = {
    0, 0, 0,
    1, 1, 1,
    2, 2, 2
};

static const uint8_t turns_table[MOVES] = {
    1, 2, 3,
    1, 2, 3,
    1, 2, 3
};

static const uint8_t inverse_move[MOVES] = {
    2, 1, 0,
    5, 4, 3,
    8, 7, 6
};

static const uint8_t mod3_table[5] = {
    0, 1, 2, 0, 1
};

/*
 * transition tables
 *
 * perm_transition[perm_rank][move]
 * ori_transition[ori_rank][move]
 */
static uint16_t perm_transition[PERMUTATIONS][MOVES];
static uint16_t ori_transition[ORIENTATIONS][MOVES];


/* ----------------------------------------------------------
 * Move implementation
 * ---------------------------------------------------------- */

static state_t quarter_turn(state_t state, uint8_t face)
{
    state_t result;

    for (uint8_t i = 0; i < CUBIES; ++i) {
        uint8_t from = source[face][i];

        result.p[i] = state.p[from];

        result.o[i] =
            mod3_table[state.o[from] + twist[face][i]];
    }

    return result;
}

static state_t apply_move(state_t state, uint8_t move)
{
    uint8_t face = face_table[move];
    uint8_t turns = turns_table[move];

    for (uint8_t t = 0; t < turns; ++t)
        state = quarter_turn(state, face);

    return state;
}


/* ----------------------------------------------------------
 * Ranking
 * ---------------------------------------------------------- */

static uint32_t permutation_rank(const state_t *state)
{
    uint32_t rank = 0;

    for (uint8_t i = 0; i < CUBIES; ++i) {
        uint8_t smaller = 0;

        for (uint8_t j = (uint8_t)(i + 1); j < CUBIES; ++j) {
            if (state->p[j] < state->p[i])
                ++smaller;
        }

        rank = rank * (CUBIES - i) + smaller;
    }

    return rank;
}

static uint32_t orientation_rank(const state_t *state)
{
    uint32_t rank = 0;

    for (uint8_t i = 0; i < 6; ++i)
        rank = rank * 3 + state->o[i];

    return rank;
}


/* ----------------------------------------------------------
 * Unranking
 * ---------------------------------------------------------- */

/*
 * Convert permutation rank 0..5039 back to p[0..6].
 */
static void unrank_permutation(uint32_t rank, state_t *state)
{
    uint8_t available[CUBIES] = {
        0, 1, 2, 3, 4, 5, 6
    };

    uint32_t p = rank;
    uint32_t factorial = 720; /* 6! */

    for (uint8_t i = 0; i < CUBIES; ++i) {
        uint8_t q = (uint8_t)(p / factorial);

        p %= factorial;

        state->p[i] = available[q];

        /*
         * Remove available[q].
         */
        for (uint8_t j = q;
             j + 1 < (uint8_t)(CUBIES - i);
             ++j) {
            available[j] = available[j + 1];
        }

        /*
         * 720 -> 120 -> 24 -> 6 -> 2 -> 1
         */
        if (i < 5)
            factorial /= (6U - i);
    }

    /*
     * Orientation does not matter for permutation transitions.
     * Keep it as a legal solved orientation.
     */
    for (uint8_t i = 0; i < CUBIES; ++i)
        state->o[i] = 0;
}


/*
 * Convert orientation rank 0..728 back to o[0..6].
 *
 * Only o[0..5] are explicitly encoded.
 * o[6] follows from:
 *
 * (o0 + o1 + ... + o6) mod 3 = 0
 */
static void unrank_orientation(uint32_t rank, state_t *state)
{
    uint32_t value = rank;
    uint8_t sum = 0;

    /*
     * Decode six base-3 digits backwards.
     */
    for (int i = 5; i >= 0; --i) {
        state->o[i] = (uint8_t)(value % 3U);
        value /= 3U;

        sum = (uint8_t)(sum + state->o[i]);
    }

    /*
     * Choose seventh orientation so total twist is 0 mod 3.
     */
    state->o[6] =
        (uint8_t)((3U - (sum % 3U)) % 3U);

    /*
     * Permutation does not matter for orientation transitions.
     */
    for (uint8_t i = 0; i < CUBIES; ++i)
        state->p[i] = i;
}


/* ----------------------------------------------------------
 * Generate transition tables
 * ---------------------------------------------------------- */

static void generate_perm_transition(void)
{
    for (uint32_t rank = 0;
         rank < PERMUTATIONS;
         ++rank) {

        state_t state;
        unrank_permutation(rank, &state);

        /*
         * Sanity check:
         * unrank(rank) should rank back to rank.
         */
        if (permutation_rank(&state) != rank) {
            fprintf(stderr,
                    "permutation unrank failed at rank %u\n",
                    (unsigned)rank);
            exit(EXIT_FAILURE);
        }

        for (uint8_t move = 0; move < MOVES; ++move) {
            state_t next = apply_move(state, move);

            uint32_t next_rank =
                permutation_rank(&next);

            if (next_rank >= PERMUTATIONS) {
                fprintf(stderr,
                        "invalid permutation rank %u\n",
                        (unsigned)next_rank);
                exit(EXIT_FAILURE);
            }

            perm_transition[rank][move] =
                (uint16_t)next_rank;
        }
    }
}


static void generate_ori_transition(void)
{
    for (uint32_t rank = 0;
         rank < ORIENTATIONS;
         ++rank) {

        state_t state;
        unrank_orientation(rank, &state);

        /*
         * Sanity check:
         * unrank(rank) should rank back to rank.
         */
        if (orientation_rank(&state) != rank) {
            fprintf(stderr,
                    "orientation unrank failed at rank %u\n",
                    (unsigned)rank);
            exit(EXIT_FAILURE);
        }

        for (uint8_t move = 0; move < MOVES; ++move) {
            state_t next = apply_move(state, move);

            uint32_t next_rank =
                orientation_rank(&next);

            if (next_rank >= ORIENTATIONS) {
                fprintf(stderr,
                        "invalid orientation rank %u\n",
                        (unsigned)next_rank);
                exit(EXIT_FAILURE);
            }

            ori_transition[rank][move] =
                (uint16_t)next_rank;
        }
    }
}


/* ----------------------------------------------------------
 * Verification
 * ---------------------------------------------------------- */

static void verify_perm_transition(void)
{
    /*
     * Applying move then inverse(move) must return
     * to the same permutation rank.
     */
    for (uint32_t rank = 0;
         rank < PERMUTATIONS;
         ++rank) {

        for (uint8_t move = 0;
             move < MOVES;
             ++move) {

            uint16_t next =
                perm_transition[rank][move];

            uint16_t back =
                perm_transition[next][inverse_move[move]];

            if (back != rank) {
                fprintf(stderr,
                        "perm transition verification failed: "
                        "rank=%u move=%u next=%u back=%u\n",
                        (unsigned)rank,
                        (unsigned)move,
                        (unsigned)next,
                        (unsigned)back);

                exit(EXIT_FAILURE);
            }
        }
    }
}


static void verify_ori_transition(void)
{
    /*
     * Same test for orientation coordinates.
     */
    for (uint32_t rank = 0;
         rank < ORIENTATIONS;
         ++rank) {

        for (uint8_t move = 0;
             move < MOVES;
             ++move) {

            uint16_t next =
                ori_transition[rank][move];

            uint16_t back =
                ori_transition[next][inverse_move[move]];

            if (back != rank) {
                fprintf(stderr,
                        "orientation transition verification failed: "
                        "rank=%u move=%u next=%u back=%u\n",
                        (unsigned)rank,
                        (unsigned)move,
                        (unsigned)next,
                        (unsigned)back);

                exit(EXIT_FAILURE);
            }
        }
    }
}


/* ----------------------------------------------------------
 * Assembly output
 * ---------------------------------------------------------- */

static void print_perm_transition(void)
{
    /*
     * .half needs 2-byte alignment for convenient RV32I lhu.
     */
    puts(".balign 2");
    puts("perm_transition:");

    for (uint32_t rank = 0;
         rank < PERMUTATIONS;
         ++rank) {

        printf("    .half ");

        for (uint8_t move = 0;
             move < MOVES;
             ++move) {

            printf("%u",
                   (unsigned)perm_transition[rank][move]);

            if (move + 1 != MOVES)
                putchar(',');
        }

        putchar('\n');
    }
}


static void print_ori_transition(void)
{
    puts("");
    puts(".balign 2");
    puts("ori_transition:");

    for (uint32_t rank = 0;
         rank < ORIENTATIONS;
         ++rank) {

        printf("    .half ");

        for (uint8_t move = 0;
             move < MOVES;
             ++move) {

            printf("%u",
                   (unsigned)ori_transition[rank][move]);

            if (move + 1 != MOVES)
                putchar(',');
        }

        putchar('\n');
    }
}


/* ----------------------------------------------------------
 * Main
 * ---------------------------------------------------------- */

int main(void)
{
    fprintf(stderr,
            "Generating permutation transition table...\n");

    generate_perm_transition();

    fprintf(stderr,
            "Generating orientation transition table...\n");

    generate_ori_transition();

    fprintf(stderr,
            "Verifying permutation transitions...\n");

    verify_perm_transition();

    fprintf(stderr,
            "Verifying orientation transitions...\n");

    verify_ori_transition();

    fprintf(stderr,
            "All transition-table tests passed.\n");

    fprintf(stderr,
            "Permutation table: %u bytes\n",
            (unsigned)sizeof perm_transition);

    fprintf(stderr,
            "Orientation table: %u bytes\n",
            (unsigned)sizeof ori_transition);

    fprintf(stderr,
            "Total transition data: %u bytes\n",
            (unsigned)(sizeof perm_transition +
                       sizeof ori_transition));

    /*
     * Assembly goes to stdout.
     */
    print_perm_transition();
    print_ori_transition();

    return 0;
}