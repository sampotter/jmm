#include <jmm/eik2mp.h>

#include <stdio.h>

int main(int argc, char const *argv[]) {
  if (argc != 3) {
    printf("usage: %s <verts_path> <faces_path>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  char const *verts_path = argv[1];
  char const *faces_path = argv[2];

  mesh22_s *mesh;
  mesh22_alloc(&mesh);
  mesh22_init_from_binary_files(mesh, verts_path, faces_path);

  eik2mp_s *eik;
  eik2mp_alloc(&eik);
  eik2mp_init(eik, mesh);
  eik2mp_solve(eik);
  eik2mp_deinit(eik);
  eik2mp_dealloc(&eik);

  mesh22_deinit(mesh);
  mesh22_dealloc(&mesh);
}
