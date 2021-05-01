#ifndef LIBMOLFILE_PLUGIN_H
#define LIBMOLFILE_PLUGIN_H
#include "vmdplugin.h"

#ifdef __cplusplus
extern "C" {
#endif

extern int molfile_dcdplugin_init(void);
extern int molfile_dcdplugin_register(void *, vmdplugin_register_cb);
extern int molfile_dcdplugin_fini(void);
extern int molfile_jsplugin_init(void);
extern int molfile_jsplugin_register(void *, vmdplugin_register_cb);
extern int molfile_jsplugin_fini(void);
extern int molfile_pdbplugin_init(void);
extern int molfile_pdbplugin_register(void *, vmdplugin_register_cb);
extern int molfile_pdbplugin_fini(void);
extern int molfile_psfplugin_init(void);
extern int molfile_psfplugin_register(void *, vmdplugin_register_cb);
extern int molfile_psfplugin_fini(void);
extern int molfile_namdbinplugin_init(void);
extern int molfile_namdbinplugin_register(void *, vmdplugin_register_cb);
extern int molfile_namdbinplugin_fini(void);

#define MOLFILE_INIT_ALL \
    molfile_dcdplugin_init(); \
    molfile_jsplugin_init(); \
    molfile_pdbplugin_init(); \
    molfile_psfplugin_init(); \
    molfile_namdbinplugin_init(); \

#define MOLFILE_REGISTER_ALL(v, cb) \
    molfile_dcdplugin_register(v, cb); \
    molfile_jsplugin_register(v, cb); \
    molfile_pdbplugin_register(v, cb); \
    molfile_psfplugin_register(v, cb); \
    molfile_namdbinplugin_register(v, cb); \

#define MOLFILE_FINI_ALL \
    molfile_dcdplugin_fini(); \
    molfile_jsplugin_fini(); \
    molfile_pdbplugin_fini(); \
    molfile_psfplugin_fini(); \
    molfile_namdbinplugin_fini(); \

#ifdef __cplusplus
}
#endif
#endif
