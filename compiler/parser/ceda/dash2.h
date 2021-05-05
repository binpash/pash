struct stackmark* Dash_init_stack (void);
void Dash_pop_stack (struct stackmark* smark);
char* Dash_alloc_stack_string (char* str);
void Dash_free_stack_string (char* str);
void Dash_dash_init (void);
void Dash_initialize_dash_errno (void);
void Dash_initialize (void);
void Dash_popfile (void);
void Dash_setinputstring (char* str);
void Dash_setinputtostdin (void);
int Dash_setinputfile (char* str, int push);
void* Dash_setvar (const char* name, const char* val);
void Dash_setalias (const char *name, const char *val);
void Dash_unalias (const char* name);
int Dash_freshfd_ge10 (int fd);

union node* Dash_parsecmd_safe (int interact);
union node* Dash_parse_next (int interactive);
