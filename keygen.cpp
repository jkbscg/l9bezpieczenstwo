#include <cstdlib>
#include <cstring>
#include <gmp.h>
#include <ctime>
#include <iostream>
#include <sys/time.h>

using namespace std;

bool is_safe_prime(mpz_t p) {
    // p = 2q + 1
    mpz_t q; mpz_init(q);
    mpz_set(q, p);
    mpz_sub_ui(q, q, 1);
    mpz_tdiv_q_ui(q, q, 2);
    bool ret = mpz_probab_prime_p(q, 24);
    mpz_clear(q);
    return ret;
}

void key_gen(mpz_t p, mpz_t q, mpz_t N, int length, gmp_randstate_t state) {
    mpz_t low, high, phi_N;
    mpz_init(low); mpz_init(high); mpz_init(phi_N);
    
    mpz_ui_pow_ui(low, 2, (length/2) - 1);

    while (mpz_cmp(low, p) > 0 || mpz_probab_prime_p(p, 24) < 1) {
        mpz_urandomb(p, state, (length/2)+1);
    }

    mpz_mul(low, low, low); 
    mpz_mul_ui(low, low, 2); 
    mpz_mul_ui(high, low, 2); 
    mpz_cdiv_q(low, low, p);
    mpz_fdiv_q(high, high, p);

    while (mpz_cmp(low, q) >= 0 || mpz_cmp(q, high) > 0 || mpz_probab_prime_p(q, 24) < 1){
        mpz_urandomb(q, state, (length/2) + 1);
    }

    mpz_mul(N, p, q);
    mpz_sub_ui(p, p, 1);
    mpz_sub_ui(q, q, 1);
    mpz_mul(phi_N, p, q);

    do {
        mpz_urandomm(p, state, phi_N);
        mpz_gcd(q, phi_N, p);
    } while (mpz_cmp_ui(q,1) != 0);

    mpz_invert(q, p, phi_N);
    mpz_clear(low); mpz_clear(high); mpz_clear(phi_N);
}

void key_gen_safe(mpz_t p, mpz_t q, mpz_t N, int length, gmp_randstate_t state) {
    mpz_t low, high, phi_N;
    mpz_init(low); mpz_init(high); mpz_init(phi_N);

    mpz_ui_pow_ui(low, 2, (length/2) - 1);

    while (mpz_cmp(low, p) > 0 || mpz_probab_prime_p(p, 24) < 1 || !is_safe_prime(p)) {
        mpz_urandomb(p, state, (length/2)+1);
    }

    mpz_mul(low, low, low); 
    mpz_mul_ui(low, low, 2); 
    mpz_mul_ui(high, low, 2); 
    mpz_cdiv_q(low, low, p);
    mpz_fdiv_q(high, high, p);

    while (mpz_cmp(low, q) >= 0 || mpz_cmp(q, high) > 0 || mpz_probab_prime_p(q, 24) < 1 || 
    !is_safe_prime(q)){
        mpz_urandomb(q, state, (length/2) + 1);
    }

    mpz_mul(N, p, q);
    // mpz_out_str(stdout, 10, p); cout << ' ';
    // mpz_out_str(stdout, 10, q); cout << ' ';
    mpz_sub_ui(p, p, 1);
    mpz_sub_ui(q, q, 1);
    mpz_mul(phi_N, p, q);

    do {
        mpz_urandomm(p, state, phi_N);
        mpz_gcd(q, phi_N, p);
    } while (mpz_cmp_ui(q,1) != 0);

    mpz_invert(q, p, phi_N);
    mpz_clear(low); mpz_clear(high); mpz_clear(phi_N);

}

// p-1 pollarda
void key_gen_bad(mpz_t p, mpz_t q, mpz_t N, int length, gmp_randstate_t state) {
    mpz_t low, high, phi_N;
    mpz_init(low); mpz_init(high); mpz_init(phi_N);

    mpz_set_si(p, 65537);


    mpz_ui_pow_ui(high, 2, length);
    mpz_cdiv_q_ui(low, high, 2);
    mpz_cdiv_q(low, low, p);
    mpz_fdiv_q(high, high, p);

    while (mpz_cmp(low, q) >= 0 || mpz_cmp(q, high) > 0 || mpz_probab_prime_p(q, 24) < 1 || 
    !is_safe_prime(q)){
        mpz_urandomb(q, state, length-16);
    }

    mpz_mul(N, p, q);
    // mpz_out_str(stdout, 10, p); cout << ' ';
    // mpz_out_str(stdout, 10, q); cout << ' ';
    mpz_sub_ui(p, p, 1);
    mpz_sub_ui(q, q, 1);
    mpz_mul(phi_N, p, q);

    do {
        mpz_urandomm(p, state, phi_N);
        mpz_gcd(q, phi_N, p);
    } while (mpz_cmp_ui(q,1) != 0);

    mpz_invert(q, p, phi_N);
    mpz_clear(low); mpz_clear(high); mpz_clear(phi_N);
}




int main(int argc, char const *argv[]){
    
    string mode = argv[1];

    int length = atoi(argv[2]);
    mpz_t e, d, N;
    mpz_init(e); mpz_init(d); mpz_init(N);

    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, time(NULL));

    struct timeval diff, startTV, endTV;

    if(!mode.compare("normal")) {
        gettimeofday(&startTV, NULL); 
        key_gen(e,d,N,length,state);
        gettimeofday(&endTV, NULL); 
    }
    else if(!mode.compare("safe")) {
        gettimeofday(&startTV, NULL); 
        key_gen_safe(e,d,N,length,state);
        gettimeofday(&endTV, NULL); 
    }
    else if(!mode.compare("bad")) {
        gettimeofday(&startTV, NULL); 
        key_gen_bad(e,d,N,length,state);
        gettimeofday(&endTV, NULL); 
    }
    else {
        cout << "Złe dane"<< endl;
    }

    timersub(&endTV, &startTV, &diff);
    unsigned long time_in_micros = 1000000 * diff.tv_sec + diff.tv_usec;
    cout<< time_in_micros<<" in microseconds"<< endl;
    
    mpz_out_str(stdout, 10, d); cout << ' ';
    mpz_out_str(stdout, 10, e); cout << ' ';
    mpz_out_str(stdout, 10, N); cout << ' ';
}