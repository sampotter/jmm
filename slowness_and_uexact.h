double slowness(char slo_fun,struct myvector x);
struct myvector gradslo(char slo_fun,struct myvector x,double s);
struct mymatrix Hslo(char slo_fun,struct myvector x,struct myvector gs,double s);
double exact_solution(char slo_fun,struct myvector x,double s);
struct myvector exact_gradient(char slo_fun,struct myvector x,double s);
void set_params(char slo_fun);
