const char *argp_program_version = "3d_wedge";
char const *argp_program_bug_address = "sfp@cims.nyu.edu";

static char doc[] = ""; // TODO: add some docs here when it makes sense

static char args_doc[] = "EXAMPLE";

static struct argp_option options[] = {
  {"verbose", 'v', 0, 0, "Produce verbose output"},
  {0}
};

struct arguments {
  char *args[2];
  bool verbose;
  char *output_file;
};

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
  struct arguments *arguments = state->input;

  switch (key) {
  case 'v':
    arguments->verbose = true;
    break;
  case 'o':
    arguments->output_file = arg;
    break;
  case ARGP_KEY_ARG:
    if (state->arg_num >= 2) // too many arguments
      argp_usage(state);
    arguments->args[state->arg_num] = arg;
    break;
  case ARGP_KEY_END:
    if (state->arg_num < 2) // too few arguments
      argp_usage(state);
    break;
  default:
    return ARGP_ERR_UNKNOWN;
  }

  return EXIT_SUCCESS;
}

static struct argp argp = {options, parse_opt, args_doc, doc};

int main(int argc, char *argv[]) {
  struct arguments arguments = {
    .verbose = false,
    .output_file = "-"
  };

  argp_parse(&argp, argc, argv, 0, 0, &arguments);


  exit(EXIT_SUCCESS);
}
