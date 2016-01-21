#ifdef EXPOSE_FUNCTION_SIGNATURES
namespace asttooling {
    void ast_tooling__testDerivable(clang::ast_matchers::MatchFinder::MatchCallback *ptr);
    core::T_mv ast_tooling__deduplicate(core::List_sp replacements);
    core::Str_sp ast_tooling__clangVersionString();
    AstVisitor_sp ast_tooling__makeAstVisitor(core::T_sp target);
};
namespace serveEvent {
    core::Integer_mv serve_event_internal__ll_serveEventWithTimeout(core::ForeignData_sp rfd, core::ForeignData_sp wfd, int maxfdp1, double seconds);
    core::Integer_mv serve_event_internal__ll_serveEventNoTimeout(core::ForeignData_sp rfd, core::ForeignData_sp wfd, int maxfdp1);
    int serve_event_internal__ll_fdset_size();
    int serve_event_internal__ll_fd_isset(int fd, core::ForeignData_sp fdset);
    void serve_event_internal__ll_fd_set(int fd, core::ForeignData_sp fdset);
    void serve_event_internal__ll_fd_zero(core::ForeignData_sp fdset);
};
namespace sockets {
    bool sockets_internal__ll_setSockoptLinger(int fd, int level, int constant, int value);
    bool sockets_internal__ll_setSockoptTimeval(int fd, int level, int constant, double value);
    bool sockets_internal__ll_setSockoptBool(int fd, int level, int constant, bool value);
    bool sockets_internal__ll_setSockoptInt(int fd, int level, int constant, int value);
    core::Integer_sp sockets_internal__ll_getSockoptLinger(int fd, int level, int constant);
    core::DoubleFloat_sp sockets_internal__ll_getSockoptTimeval(int fd, int level, int constant);
    core::T_sp sockets_internal__ll_getSockoptBool(int fd, int level, int constant);
    core::Integer_sp sockets_internal__ll_getSockoptInt(int fd, int level, int constant);
    core::Str_sp sockets_internal__ll_strerror_errno();
    core::Str_sp sockets_internal__ll_strerror(int num);
    void sockets_internal__ll_autoCloseTwoWayStream(core::Stream_sp stream);
    core::T_sp sockets_internal__ll_makeStreamFromFd(const string &name,
                            int fd,
                            int streamMode,
                            core::T_sp elementType,
                            core::T_sp externalFormat);
    int sockets_internal__ll_dup(int fd);
    int sockets_internal__ll_setfNonBlockingMode(int fd, int nblock);
    int sockets_internal__ll_nonBlockingMode(int fd);
    core::Str_sp sockets_internal__socketPeername_localSocket(int fd);
    int sockets_internal__ll_socketConnect_localSocket(int fd, int family, const string &path);
    core::T_mv sockets_internal__ll_socketAccept_localSocket(int socketFileDescriptor);
    int sockets_internal__ll_socketBind_localSocket(int fd, const string &name, int family);
    core::Integer_sp sockets_internal__socketSendNoAddress(int fb,
                                  core::Vector_sp vbuffer,
                                  int length,
                                  bool oob,
                                  bool eor,
                                  bool dontroute,
                                  bool dontwait,
                                  bool nosignal,
                                  bool confirm
                                  );
    core::Integer_sp sockets_internal__ll_socketSendAddress(int fd,
                                   core::Vector_sp vbuffer,
                                   int length,
                                   int secondAddress,
                                   int ip0,
                                   int ip1,
                                   int ip2,
                                   int ip3,
                                   bool oob,
                                   bool eor,
                                   bool dontroute,
                                   bool dontwait,
                                   bool nosignal,
                                   bool confirm);
    int sockets_internal__ll_socketName(int fd, core::Vector_sp vector);
    int sockets_internal__ll_socketPeername_inetSocket(int fd, core::Vector_sp vector);
    int sockets_internal__ll_socketConnect_inetSocket(int port, int ip0, int ip1, int ip2, int ip3, int socket_file_descriptor);
    core::T_mv sockets_internal__ll_socketAccept_inetSocket(int sfd);
    int sockets_internal__ll_socketBind_inetSocket(int port, int ip0, int ip1, int ip2, int ip3, int socketFileDescriptor);
    int sockets_internal__ll_getProtocolByName(const string &name);
    core::T_mv sockets_internal__ll_socketReceive(int fd,
                         core::T_sp buffer,
                         int length,
                         bool oob,
                         bool peek,
                         bool waitall);
    core::T_sp sockets_internal__ll_getHostByAddress(core::Vector_sp address,
                            core::T_sp tHostEnt,
                            core::T_sp setf_host_ent_name,
                            core::T_sp setf_host_ent_aliases,
                            core::T_sp setf_host_ent_address_type,
                            core::T_sp setf_host_ent_addresses);
    core::T_sp sockets_internal__ll_getHostByName(const string &hostName,
                                             core::T_sp tHostEnt,
                                             core::Function_sp setf_host_ent_name,
                                             core::Function_sp setf_host_ent_aliases,
                                             core::Function_sp setf_host_ent_address_type,
                                             core::Function_sp setf_host_ent_addresses);
    core::Str_sp sockets_internal__ll_getNameServiceErrorMessage(int num);
    int sockets_internal__ll_get_name_service_h_errno();
    int sockets_internal__ll_socket_errno();
    core::T_mv sockets_internal__funcName();
    int sockets_internal__ff_close(int filedes);
    int sockets_internal__ff_listen(int sockfd, int backlog);
    int sockets_internal__ff_socket(int domain, int type, int protocol);
};
namespace gctools {
    core::Fixnum gctools__next_header_kind();
    void gctools__debug_allocations(core::T_sp debugOn);
    void gctools__cleanup();
    void gctools__garbage_collect();
    core::T_sp gctools__stack_depth();
    core::T_mv cl__room(core::T_sp x, core::Fixnum_sp marker, core::T_sp tmsg);
    core::Symbol_sp gctools__alloc_pattern_end();
    void gctools__alloc_pattern_begin(core::Symbol_sp pattern);
    Fixnum gctools__gc_marker(core::Fixnum_sp marker);
    void gctools__monitor_allocations(bool on, core::Fixnum_sp backtraceStart, core::Fixnum_sp backtraceCount, core::Fixnum_sp backtraceDepth);
    core::T_mv gctools__gc_info(core::T_sp x, core::Fixnum_sp marker);
    core::T_mv core__hardwired_kinds();
    bool core__inherits_from_instance(core::T_sp obj);
    Fixnum core__header_kind(core::T_sp obj);
    core::T_mv gctools__bytes_allocated();
    void gctools__deallocate_unmanaged_instance(core::T_sp obj);
    core::T_sp gctools__bootstrap_kind_p(const string &name);
    core::Cons_sp gctools__bootstrap_kind_symbols();
    int gctools__max_bootstrap_kinds();
};
namespace telemetry {
    size_t core__telemetry_count();
    void core__telemetry_dump(core::T_sp begin, core::T_sp end);
    void core__telemetry_labels();
    void core__telemetry_follow(core::T_sp address);
    void core__telemetry_search_labels(core::List_sp labels);
    void core__telemetry_search(core::List_sp addresses);
    void core__telemetry_open(core::T_sp pathname);
};
namespace llvmo {
    void llvm_sys__viewCFG(core::T_sp funcDes, core::T_sp only);
    void llvm_sys__disassembleSTAR(core::Function_sp cf);
    llvmo::GlobalVariable_sp llvm_sys__getOrCreateExternalGlobal(llvmo::Module_sp module, const string &name, llvmo::Type_sp data_type);
    void llvm_sys__throwIfMismatchedStructureSizes(core::Fixnum_sp tspSize, core::Fixnum_sp tmvSize, gc::Nilable<core::Fixnum_sp> givenIhfSize);
    core::T_sp llvm_sys__cxxDataStructuresInfo();
    core::Str_sp llvm_sys__mangleSymbolName(core::Str_sp name);
    core::T_mv TargetRegistryLookupTarget(const std::string &ArchName, Triple_sp triple);
    void finalizeClosure(ExecutionEngine_sp oengine, core::Function_sp func);
    core::Function_sp finalizeEngineAndRegisterWithGcAndGetCompiledFunction(ExecutionEngine_sp oengine, core::T_sp functionName, Function_sp fn, core::T_sp activationFrameEnvironment, core::Str_sp globalRunTimeValueName, core::T_sp fileName, size_t filePos, int linenumber, core::T_sp lambdaList);
    Function_sp llvm_sys__FunctionCreate(FunctionType_sp tysp, llvm::GlobalValue::LinkageTypes linkage, core::Str_sp nsp, Module_sp modulesp);
    void PassManagerBuilderSetfSizeLevel(PassManagerBuilder_sp pmb, int level);
    void PassManagerBuilderSetfOptLevel(PassManagerBuilder_sp pmb, int optLevel);
    void PassManagerBuilderSetfInliner(PassManagerBuilder_sp pmb, llvm::Pass *inliner);
    core::List_sp llvm_sys__module_get_function_list(Module_sp module);
    bool llvm_sys__valid(core::T_sp value);
    Value_sp llvm_sys__makeStringGlobal(Module_sp module, core::Str_sp svalue);
    bool llvm_sys__valuep(core::T_sp arg);
    Module_sp llvm_sys__parseBitcodeFile(core::Str_sp filename, LLVMContext_sp context);
    void llvm_sys__writeBitcodeToFile(Module_sp module, core::Str_sp pathname);
    core::T_mv llvm_sys__verifyFunction(Function_sp function);
    core::T_mv llvm_sys__verifyModule(Module_sp module, core::Symbol_sp action);
    void llvm_sys__writeIrToFile(Module_sp module, core::Str_sp path);
    core::T_mv llvm_sys__link_in_module(Linker_sp linker, Module_sp module);
    bool llvm_sys__llvm_value_p(core::T_sp o);
    void compiler__setAssociatedFuncs(core::CompiledFunction_sp cf, core::List_sp associatedFuncs);
    uint llvm_sys__dwTag(core::Symbol_sp tagsym, uint debugVersion);
    void llvm_sys__addGlobalBootFunctionsSizePass(llvmo::PassManager_sp passManager);
};
namespace cffi {
    Pointer_sp cffi_sys__foreign_alloc(core::Integer_sp size);
    core::Fixnum_sp cffi_sys__PERCENTforeign_type_size(core::Symbol_sp atype);
    core::Fixnum_sp cffi_sys__PERCENTforeign_type_alignment(core::Symbol_sp atype);
    Pointer_sp cffi_sys__foreign_symbol_pointer(core::Str_sp name);
    Pointer_sp cffi_sys__PERCENTload_foreign_library(core::Str_sp name);
};
namespace core {
    T_sp cl__get_internal_run_time();
    T_sp cl__get_internal_real_time();
    T_mv cl__ensure_directories_exist(T_sp pathspec);
    void core__load_time_values_dump_symbols(T_sp nameOrLtv, T_sp indices);
    void core__load_time_values_dump_values(T_sp nameOrLtv, T_sp indices);
    void core__load_time_values_ids();
    Symbol_sp core__lookup_load_time_symbol(const string &name, int idx);
    T_sp core__lookup_load_time_value(const string &name, int idx);
    LoadTimeValues_mv core__load_time_value_array(const string &name, int dataSize, int symbolSize);
    bool core__set_run_time_values_vector(const string &name);
    Fixnum_sp cl__vector_push_extend(T_sp newElement, Vector_sp vec, int extension);
    T_sp cl__vector_push(T_sp newElement, Vector_sp vec);
    T_sp core__adjust_vector(T_sp array, int new_dimensions, T_sp initial_element, List_sp initial_contents);
    Vector_sp core__make_vector(T_sp element_type,
                           int dimension,
                           bool adjustable,
                           T_sp fill_pointer,
                           T_sp displaced_to,
                           T_sp displaced_index_offset,
                           T_sp initial_element,
                           T_sp initialContents);
    Vector_sp cl__vector(List_sp args);
    T_sp core__setf_subseq(T_sp sequence, int start, Fixnum_sp end, T_sp subseq);
    T_sp core__elt_set(T_sp sequence, int index, T_sp val);
    T_sp core__setf_elt(T_sp sequence, int index, T_sp value);
    T_sp cl__copy_seq(T_sp seq);
    T_sp cl__subseq(T_sp seq, int start, T_sp end);
    T_sp cl__nreverse(T_sp seq);
    T_sp cl__reverse(T_sp seq);
    T_mv cl__make_sequence(T_sp type, Fixnum_sp size, T_sp initial_element, T_sp iesp);
    T_sp cl__elt(T_sp sequence, int index);
    uint cl__length(T_sp arg);
    core::Package_sp coerce_to_package(core::T_sp obj);
    core::Path_sp core__path_designator(core::T_sp obj);
    T_mv cl__parse_integer(Str_sp str, Fixnum start, T_sp end, uint radix, T_sp junkAllowed);
    claspChar core__schar_set(Str_sp str, int idx, claspChar c);
    claspChar core__char_set(Str_sp str, int idx, claspChar c);
    claspChar cl__schar(Str_sp str, int idx);
    T_mv cl__make_string(Fixnum_sp size, T_sp initial_element, T_sp element_type);
    T_mv cl__string_not_lessp(T_sp strdes1, T_sp strdes2, Fixnum_sp start1, T_sp end1, Fixnum_sp start2, T_sp end2);
    T_mv cl__string_not_greaterp(T_sp strdes1, T_sp strdes2, Fixnum_sp start1, T_sp end1, Fixnum_sp start2, T_sp end2);
    T_mv cl__string_greaterp(T_sp strdes1, T_sp strdes2, Fixnum_sp start1, T_sp end1, Fixnum_sp start2, T_sp end2);
    T_mv cl__string_lessp(T_sp strdes1, T_sp strdes2, Fixnum_sp start1, T_sp end1, Fixnum_sp start2, T_sp end2);
    T_mv cl__string_not_equal(T_sp strdes1, T_sp strdes2, Fixnum_sp start1, T_sp end1, Fixnum_sp start2, T_sp end2);
    T_sp cl__string_equal(T_sp strdes1, T_sp strdes2, Fixnum_sp start1, T_sp end1, Fixnum_sp start2, T_sp end2);
    T_mv cl__string_GE_(T_sp strdes1, T_sp strdes2, Fixnum_sp start1, T_sp end1, Fixnum_sp start2, T_sp end2);
    T_mv cl__string_LE_(T_sp strdes1, T_sp strdes2, Fixnum_sp start1, T_sp end1, Fixnum_sp start2, T_sp end2);
    T_mv cl__string_GT_(T_sp strdes1, T_sp strdes2, Fixnum_sp start1, T_sp end1, Fixnum_sp start2, T_sp end2);
    T_mv cl__string_LT_(T_sp strdes1, T_sp strdes2, Fixnum_sp start1, T_sp end1, Fixnum_sp start2, T_sp end2);
    T_mv cl__string_NE_(T_sp strdes1, T_sp strdes2, Fixnum_sp start1, T_sp end1, Fixnum_sp start2, T_sp end2);
    T_sp cl__string_EQ_(T_sp strdes1, T_sp strdes2, Fixnum_sp start1, T_sp end1, Fixnum_sp start2, T_sp end2);
    Str_sp cl__string(T_sp arg);
    Str_sp cl__nstring_downcase(T_sp arg);
    Str_sp cl__nstring_upcase(T_sp arg);
    Str_sp cl__string_downcase(T_sp arg);
    claspChar cl__char(T_sp ostr, int idx);
    Str_sp cl__string_upcase(T_sp arg);
    Str_sp cl__string_right_trim(T_sp charbag, T_sp str);
    Str_sp cl__string_left_trim(T_sp charbag, T_sp str);
    Str_sp cl__string_trim(T_sp charbag, T_sp str);
    Character_sp cl__code_char(Integer_sp ich);
    Fixnum_sp cl__char_int(Character_sp och);
    Fixnum_sp cl__char_code(Character_sp och);
    Str_sp cl__char_name(Character_sp och);
    T_mv cl__name_char(Str_sp sname);
    T_sp cl__digit_char_p(Character_sp c, Fixnum_sp radix);
    bool cl__alpha_char_p(Character_sp ch);
    bool cl__standard_char_p(Character_sp ch);
    bool cl__char_equal(List_sp args);
    T_mv cl__char_not_equal(List_sp args);
    T_sp cl__char_EQ_(List_sp args);
    T_sp cl__char_NE_(List_sp args);
    T_sp cl__char_not_lessp(List_sp args);
    T_mv cl__char_not_greaterp(List_sp args);
    T_mv cl__char_greaterp(List_sp args);
    T_mv cl__char_lessp(List_sp args);
    T_mv cl__char_GE_(List_sp args);
    T_sp cl__char_LE_(List_sp args);
    T_sp cl__char_GT_(List_sp args);
    T_sp cl__char_LT_(List_sp args);
    char cl__char_downcase(Character_sp ch);
    char cl__char_upcase(Character_sp ch);
    bool cl__alphanumericp(Character_sp ch);
    bool cl__both_case_p(Character_sp c);
    bool cl__upper_case_p(Character_sp c);
    bool cl__lower_case_p(Character_sp c);
    bool core__structure_subtypep(T_sp x, Symbol_sp y);
    bool core__structurep(T_sp arg);
    T_sp core__structure_set(T_sp obj, Symbol_sp type, int idx, T_sp val);
    T_sp core__structure_ref(T_sp obj, Symbol_sp type, int idx);
    T_sp cl__copy_structure(T_sp arg);
    T_sp core__make_structure(T_sp type, List_sp slot_values);
    void core__ensure_single_dispatch_method(Symbol_sp gfname, Class_sp receiver_class, LambdaListHandler_sp lambda_list_handler, List_sp declares, gc::Nilable<Str_sp> docstring, Function_sp body);
    T_sp core__ensure_single_dispatch_generic_function(Symbol_sp gfname, LambdaListHandler_sp llhandler);
    T_mv core__nread(T_sp sin, T_sp eof_error_p, T_sp eof_value);
    void core__ensure_documentation(T_sp sub_key, Symbol_sp symbol, Str_sp value);
    T_mv ext__annotate(T_sp object, T_sp key, T_sp sub_key, Str_sp value);
    T_sp core__rem_record_field(List_sp record, T_sp key, T_sp sub_key);
    T_sp core__set_record_field(List_sp record, T_sp key, T_sp sub_key, Str_sp value);
    T_sp core__record_field(List_sp record, T_sp key, T_sp sub_key);
    T_sp core__record_cons(List_sp record, T_sp key, T_sp sub_key);
    T_sp core__backquote_remove_tokens(T_sp x);
    List_sp core__backquote_attach_conses(T_sp items, T_sp result);
    T_sp core__backquote_attach_append(T_sp op, T_sp item, T_sp result);
    T_sp core__backquote_null_or_quoted(T_sp x);
    T_sp core__backquote_simplify_args(T_sp x);
    T_sp core__backquote_simplify(T_sp x);
    T_sp core__backquote_maptree(Function_sp op, T_sp x);
    bool core__backquote_frob(T_sp x);
    bool core__backquote_splicing_frob(T_sp x);
    T_sp core__backquote_bracket(T_sp x);
    T_sp core__backquote_process(T_sp ox);
    T_mv core__backquote_completely_process(T_sp x);
    T_sp core__backquote_append(List_sp lists);
    T_sp cl__format(T_sp destination, T_sp control, List_sp args);
    T_sp core__bformat(T_sp destination, const string &control, List_sp args);
    T_sp core__rem_sysprop(T_sp key, T_sp area);
    T_mv core__get_sysprop(T_sp key, T_sp area);
    T_sp core__put_sysprop(T_sp key, T_sp area, T_sp value);
    Iterator_sp core__iterator_step(Iterator_sp it);
    T_sp core__iterator_unsafe_element(Iterator_sp it);
    string core__condition_message(T_sp condition);
    T_sp cl__make_condition(T_sp type, List_sp slot_initializations);
    Symbol_mv cl__make_symbol(Str_sp name);
    T_sp core__symbol_value_address(const Symbol_sp arg);
    T_sp cl__symbol_value(const Symbol_sp arg);
    Str_sp cl__symbol_name(Symbol_sp arg);
    Function_sp cl__symbol_function(Symbol_sp sym);
    T_sp cl__symbol_package(Symbol_sp arg);
    bool cl__boundp(Symbol_sp arg);
    T_sp core__putprop(Symbol_sp sym, T_sp val, T_sp indicator);
    void core__setf_symbol_plist(Symbol_sp sym, List_sp plist);
    T_sp cl__get(Symbol_sp sym, T_sp indicator, T_sp defval);
    List_sp cl__symbol_plist(Symbol_sp sym);
    string core__alist_asString(List_sp alist);
    List_sp core__alist_get(List_sp alist, T_sp key);
    List_sp core__alist_push(List_sp alist, T_sp key, T_sp val);
    List_sp core__alist_erase(List_sp alist, T_sp key);
    T_sp core__append2(List_sp x, List_sp y);
    List_sp cl__make_list(Integer_sp osize, T_sp initial_element);
    Cons_sp cl__rplacd(Cons_sp c, T_sp o);
    Cons_sp cl__rplaca(Cons_sp c, T_sp o);
    Cons_sp cl__cons(T_sp obj1, T_sp obj2);
    T_mv core__rem_f(List_sp plist, Symbol_sp indicator);
    T_sp cl__getf(List_sp plist, T_sp indicator, T_sp default_value);
    List_sp core__put_f(List_sp place, T_sp value, T_sp indicator);
    T_sp core__aset(Array_sp array, List_sp indices_value);
    void core__copy_subarray(Array_sp out, Fixnum_sp outStart, Array_sp in, Fixnum_sp inStart, Fixnum_sp len);
    T_mv cl__upgraded_array_element_type(T_sp type);
    T_mv cl__array_displacement(T_sp array);
    int cl__fill_pointer(Vector_sp vec);
    T_sp core__calculate_missing_common_lisp_symbols();
    T_sp cl__package_name(T_sp pkgDesig);
    List_sp cl__package_used_by_list(T_sp pkgDesig);
    T_sp cl__package_use_list(T_sp package_designator);
    T_mv cl__gentemp(T_sp prefix, T_sp package_designator);
    T_mv cl__shadowing_import(T_sp symbol_names_desig, T_sp package_desig);
    T_mv cl__shadow(List_sp symbol_names_desig, T_sp package_desig);
    T_mv cl__import(T_sp symbols_desig, T_sp package_desig);
    T_sp cl__package_shadowing_symbols(T_sp package_desig);
    T_sp cl__delete_package(T_sp pobj);
    T_sp cl__unuse_package(T_sp packages_to_unuse_desig, T_sp package_desig);
    T_sp cl__use_package(T_sp packages_to_use_desig, T_sp package_desig);
    T_sp cl__list_all_packages();
    T_mv cl__make_package(T_sp package_name_desig, List_sp nick_names, List_sp use_packages);
    T_mv cl__find_symbol(const string &symbolname, T_sp packageDesig);
    bool cl__unintern(Symbol_sp sym, T_sp packageDesig);
    T_sp cl__package_nicknames(T_sp pkg);
    Package_sp cl__rename_package(T_sp pkg, T_sp newNameDesig, T_sp nickNameDesigs);
    T_sp core__base_string_concatenate(T_sp args);
    T_sp core__search_string(Str_sp str1, Fixnum_sp start1, T_sp end1, Str_sp str2, Fixnum_sp start2, T_sp end2);
    bool core__logical_pathname_p(T_sp obj);
    bool cl__pathnamep(T_sp obj);
    bool core__proper_list_p(T_sp arg);
    bool core__generic_function_p(T_sp o);
    bool cl__compiled_function_p(T_sp o);
    bool cl__simple_vector_p(T_sp o);
    bool cl__simple_bit_vector_p(T_sp o);
    bool core__external_object_p(T_sp obj);
    bool core__activation_frame_p(T_sp obj);
    bool core__single_dispatch_generic_function_p(T_sp obj);
    bool core__arrayp(T_sp obj);
    bool core__structure_object_p(T_sp obj);
    bool cl__readtablep(T_sp obj);
    bool cl__hash_table_p(T_sp obj);
    bool cl__bit_vector_p(T_sp obj);
    bool core__path_p(T_sp obj);
    bool cl__keywordp(T_sp obj);
    bool cl__integerp(T_sp obj);
    bool cl__vectorp(T_sp obj);
    bool cl__floatp(T_sp obj);
    bool cl__realp(T_sp obj);
    bool core__single_float_p(T_sp obj);
    bool core__short_float_p(T_sp obj);
    bool core__long_float_p(T_sp obj);
    bool cl__rationalp(T_sp obj);
    bool cl__random_state_p(T_sp obj);
    bool core__ratio_p(T_sp obj);
    bool cl__complexp(T_sp obj);
    bool cl__numberp(T_sp obj);
    bool core__lambda_list_handler_p(T_sp obj);
    bool core__cxx_instance_p(T_sp obj);
    bool core__built_in_class_p(T_sp obj);
    bool clos__classp(T_sp obj);
    bool cl__packagep(T_sp obj);
    bool core__interpreted_function_p(T_sp arg);
    bool cl__functionp(T_sp obj);
    bool core__double_float_p(T_sp obj);
    bool core__simple_string_p(T_sp obj);
    bool cl__stringp(T_sp obj);
    bool core__bignump(T_sp obj);
    bool core__base_char_p(T_sp arg);
    T_sp core__test_val(T_sp v);
    void core__is_array(T_sp arg);
    void core__is_string(T_sp arg);
    bool core__subclassp(T_sp low, T_sp high);
    T_sp core__allocate_raw_class(T_sp orig, Class_sp metaClass, int slots, T_sp className);
    void core__inherit_default_allocator(Class_sp cl, T_sp directSuperclasses);
    void cl__reader_error(const string &sourceFile, uint lineno, Symbol_sp function,
                    Str_sp fmt, List_sp fmtargs, T_sp stream);
    void core__wrong_index(const string &sourceFile, int lineno, Symbol_sp function, T_sp array, int which, T_sp index, int nonincl_limit);
    void core__wrong_type_nth_arg(const string &sourceFile, int lineno, Symbol_sp function, int narg, T_sp value, T_sp type);
    void core__wrong_type_argument(const string &sourceFile, int lineno, Symbol_sp function, T_sp value, T_sp type);
    T_sp core__signal_simple_error(T_sp baseCondition, T_sp continueMessage, T_sp formatControl, T_sp formatArgs, T_sp args);
    Pathname_sp cl__translate_logical_pathname(T_sp tsource);
    Pathname_sp cl__translate_pathname(T_sp tsource, T_sp tfrom, T_sp tto, T_sp scase);
    T_sp core__pathname_translations(T_sp host, T_sp hostp, T_sp set);
    bool cl__pathname_match_p(T_sp tpath, T_sp tmask);
    Str_sp cl__enough_namestring(T_sp tpath, T_sp tdefaults);
    Str_sp cl__host_namestring(T_sp tpname);
    Str_sp cl__directory_namestring(T_sp tpname);
    Str_sp cl__file_namestring(T_sp tpname);
    T_sp cl__pathname_version(T_sp tpname);
    T_sp cl__pathname_type(T_sp tpname, Symbol_sp scase);
    T_sp cl__pathname_name(T_sp tpname, Symbol_sp scase);
    T_sp cl__pathname_directory(T_sp tpname, Symbol_sp scase);
    T_sp cl__pathname_device(T_sp tpname, Symbol_sp scase);
    T_sp cl__pathname_host(T_sp tpname, Symbol_sp scase);
    Pathname_sp cl__make_pathname(T_sp host, bool hostp, T_sp device, bool devicep, T_sp directory, bool directoryp, T_sp name, bool namep, T_sp type, bool typep, T_sp version, bool versionp, T_sp scase, T_sp odefaults);
    T_mv cl__parse_namestring(T_sp thing, T_sp host, T_sp tdefaults, Fixnum_sp start, T_sp end, bool junkAllowed);
    T_sp cl__namestring(T_sp x);
    Str_sp core__coerce_to_filename(T_sp pathname_orig);
    Pathname_sp core__coerce_to_physical_pathname(T_sp x);
    Pathname_sp core__coerce_to_file_pathname(T_sp tpathname);
    bool cl__wild_pathname_p(T_sp tpathname, T_sp component);
    Pathname_sp cl__merge_pathnames(T_sp path, T_sp defaults, T_sp defaultVersion);
    Pathname_sp cl__pathname(T_sp x);
    Pathname_sp core__safe_default_pathname_defaults_host_only(void);
    Pathname_sp core__safe_default_pathname_defaults(void);
    T_mv core__float_to_digits(T_sp tdigits, Float_sp number, gc::Nilable<Real_sp> position,
                          T_sp relativep);
    T_mv core__sharp_vertical_bar(T_sp sin, Character_sp ch, T_sp num);
    T_mv core__sharp_minus(T_sp sin, Character_sp ch, T_sp num);
    T_mv core__sharp_plus(T_sp sin, Character_sp ch, T_sp num);
    T_sp core__reader_feature_p(T_sp feature_test);
    T_mv core__sharp_p(T_sp sin, Character_sp ch, T_sp num);
    T_mv core__sharp_s(T_sp sin, Character_sp ch, T_sp num);
    T_mv core__sharp_a(T_sp sin, Character_sp ch, T_sp num);
    T_mv core__sharp_c(T_sp sin, Character_sp ch, T_sp num);
    T_mv core__sharp_x(T_sp sin, Character_sp ch, T_sp num);
    T_mv core__sharp_o(T_sp sin, Character_sp ch, gc::Nilable<Fixnum_sp> num);
    T_mv core__sharp_b(T_sp sin, Character_sp ch, gc::Nilable<Fixnum_sp> num);
    T_mv core__sharp_r(T_sp sin, Character_sp ch, gc::Nilable<Fixnum_sp> nradix);
    T_mv core__sharp_colon(T_sp sin, Character_sp ch, T_sp num);
    T_mv core__sharp_asterisk(T_sp sin, Character_sp ch, T_sp num);
    T_mv core__sharp_left_parenthesis(T_sp sin, Character_sp ch, T_sp tnum);
    T_sp core__sharp_single_quote(T_sp sin, Character_sp ch, T_sp num);
    T_sp core__sharp_dot(T_sp sin, Character_sp ch, T_sp num);
    T_mv core__sharp_backslash(T_sp sin, Character_sp ch, T_sp num);
    T_mv core__dispatch_macro_character(T_sp sin, Character_sp ch);
    T_mv core__reader_skip_semicolon_comment(T_sp sin, Character_sp ch);
    T_sp core__reader_quote(T_sp sin, Character_sp ch);
    T_mv core__reader_error_unmatched_close_parenthesis(T_sp sin, Character_sp ch);
    T_sp core__reader_list_allow_consing_dot(T_sp sin, Character_sp ch);
    T_sp core__reader_comma_form(T_sp sin, Character_sp ch);
    T_mv core__reader_backquoted_expression(T_sp sin, Character_sp ch);
    T_mv core__reader_double_quote_string(T_sp stream, Character_sp ch);
    T_mv cl__set_macro_character(Character_sp ch, T_sp func_desig, T_sp non_terminating_p, ReadTable_sp readtable);
    T_mv cl__get_dispatch_macro_character(Character_sp dispChar, Character_sp subChar, ReadTable_sp readtable);
    T_mv cl__set_dispatch_macro_character(Character_sp dispChar, Character_sp subChar, T_sp newFunctionDesig, ReadTable_sp readtable);
    void core__readtable_case_set(ReadTable_sp readTable, T_sp mode);
    T_sp cl__readtable_case(ReadTable_sp readTable);
    T_sp cl__copy_readtable(gc::Nilable<ReadTable_sp> fromReadTable, gc::Nilable<ReadTable_sp> toReadTable);
    T_mv cl__get_macro_character(Character_sp chr, T_sp readtable);
    T_sp cl__make_dispatch_macro_character(Character_sp ch, T_sp nonTerminatingP, ReadTable_sp readtable);
    T_sp cl__set_syntax_from_char(Character_sp toChar, Character_sp fromChar, ReadTable_sp toReadTable, gc::Nilable<ReadTable_sp> fromReadTable, T_sp fromReadTableP);
    BranchSNode_sp core__make_branch_snode();
    T_sp core__pointer_address(T_sp ptr);
    void core__pointer_delete(T_sp ptr);
    Pointer_sp core__pointer_release(T_sp ptr);
    Regex_sp core__make_regex(const string &str);
    T_sp core__write_ugly_object(T_sp obj, T_sp ostrm);
    core::Bignum_mv core__test_profile_timer(uint delay);
    core::Bignum_mv core__clock_gettime_nanoseconds();
    T_sp core__deep_copy(T_sp obj);
    bool core__sl_boundp(T_sp obj);
    T_sp core__decode(T_sp obj, core::List_sp arg);
    core::List_sp core__encode(T_sp arg);
    bool core__is_nil(T_sp arg);
    T_sp core__instancep(T_sp obj);
    T_sp core__instance_ref(T_sp obj, int idx);
    T_sp core__instance_set(T_sp obj, int idx, T_sp val);
    T_sp core__instance_sig_set(T_sp arg);
    T_sp core__instance_sig(T_sp obj);
    string core__class_name_as_string(T_sp arg);
    Class_sp core__instance_class(T_sp arg);
    T_sp core__implementation_class(T_sp arg);
    T_sp cl__copy_tree(T_sp arg);
    void core__low_level_describe(T_sp obj);
    T_sp core__print_cxx_object(T_sp obj, T_sp stream);
    bool core__fieldsp(T_sp obj);
    T_sp core__make_cxx_object(T_sp class_or_name, T_sp args);
    List_sp core__all_source_files();
    T_mv cl__not(T_sp x);
    T_mv core__repr(T_sp obj);
    T_mv core__is_sub_class_of(Class_sp tag, Class_sp mc);
    T_mv core__is_assignable_to(T_sp tag, Class_sp mc);
    void cl__cerror(T_sp cformat, T_sp eformat, List_sp arguments);
    void cl__error(T_sp datum, List_sp initializers);
    void core__invoke_internal_debugger_from_gdb();
    HashTable_sp core__single_dispatch_generic_function_table();
    void core__invoke_internal_debugger(T_sp condition);
    T_mv core__universal_error_handler(T_sp continueString, T_sp datum, List_sp initializers);
    T_mv cl__intern(Str_sp symbol_name, T_sp package_desig);
    void core__export_to_python(T_sp symbolsDesig);
    T_sp cl__export(T_sp symDes, T_sp packageDes);
    void core__debug_log_off();
    void core__debug_log_on();
    T_mv core__is_top_level_script();
    T_mv core__set_current_working_directory(Str_sp dir);
    Path_mv core__database_dir();
    Path_mv core__lisp_code_path(T_sp relPathDesignator);
    Path_mv core__library_path(T_sp relPathDesignator);
    Path_mv core__script_dir();
    T_mv core__source_line_column();
    T_mv core__source_file_name();
    T_sp cl__sort(List_sp sequence, T_sp predicate);
    List_sp core__sorted(List_sp unsorted);
    T_sp cl__apropos(Str_sp string_desig, T_sp package_desig);
    T_mv cl__macroexpand(T_sp form, T_sp env);
    T_mv cl__macroexpand_1(T_sp form, T_sp env);
    T_mv core__mpi_size();
    T_mv core__mpi_rank();
    T_mv core__mpi_enabled();
    void core__select_package(T_sp package_designator);
    T_sp cl__find_package(T_sp name_desig);
    T_mv core__find_file_in_lisp_path(Str_sp partialPath);
    Class_mv core__setf_find_class(T_sp newValue, Symbol_sp name, bool errorp, T_sp env);
    Class_mv cl__find_class(Symbol_sp symbol, bool errorp, T_sp env);
    T_mv core__getline(Str_sp prompt);
    List_sp core__member1(T_sp item, List_sp list, T_sp test, T_sp test_not, T_sp key);
    List_sp core__member_test(T_sp item, List_sp list, T_sp key, T_sp test, T_sp test_not);
    List_sp cl__member(T_sp item, T_sp tlist, T_sp key, T_sp test, T_sp test_not);
    List_sp cl__assoc(T_sp item, List_sp alist, T_sp key, T_sp test, T_sp test_not);
    List_sp cl__acons(T_sp key, T_sp val, List_sp alist);
    void core__quit(int exitValue);
    void core__exit(int exitValue);
    void core__setup_stack_monitor(T_sp warnSize, T_sp sampleSize);
    T_mv core__stack_limit();
    void core__stack_monitor();
    uint core__stack_used();
    void core__low_level_repl();
    T_sp core__list_of_all_special_operators();
    T_sp cl__boole(T_sp op, T_sp arg1, T_sp arg2);
    T_sp cl__logbitp(Integer_sp p, Integer_sp x);
    T_sp core__bit_array_op(T_sp o, T_sp tx, T_sp ty, T_sp tr);
    T_sp cl__write_sequence(T_sp seq, T_sp stream, Fixnum_sp fstart, T_sp tend);
    T_sp core__file_column(T_sp strm);
    void cl__unread_char(Character_sp ch, T_sp dstrm);
    void cl__finish_output(T_sp ostrm);
    void cl__force_output(T_sp ostrm);
    bool cl__listen(T_sp strm);
    void cl__clear_output(T_sp dstrm);
    void cl__clear_input(T_sp dstrm);
    Character_sp cl__write_char(Character_sp chr, T_sp stream);
    Integer_sp cl__write_byte(Integer_sp byte, T_sp stream);
    String_sp cl__write_line(Str_sp str, T_sp stream, Fixnum_sp start, T_sp end);
    Str_sp cl__write_string(Str_sp str, T_sp stream, int start, T_sp end);
    bool cl__fresh_line(T_sp outputStreamDesig);
    void cl__terpri(T_sp outputStreamDesig);
    T_mv cl__read_line(T_sp sin, T_sp eof_error_p, T_sp eof_value, T_sp recursive_p);
    T_mv cl__read_from_string(Str_sp content, T_sp eof_error_p, T_sp eof_value, Fixnum_sp start, T_sp end, T_sp preserve_whitespace);
    T_sp cl__read_char_no_hang(T_sp strm, T_sp eof_error_p, T_sp eof_value, T_sp recursive_p);
    T_sp cl__read_char(T_sp strm, T_sp eof_error_p, T_sp eof_value, T_sp recursive_p);
    T_sp cl__peek_char(T_sp peek_type, T_sp strm, T_sp eof_errorp, T_sp eof_value, T_sp recursive_p);
    T_sp cl__read_byte(T_sp strm, T_sp eof_error_p, T_sp eof_value);
    SourcePosInfo_sp core__input_stream_source_pos_info(T_sp strm);
    int core__stream_column(T_sp tstream);
    int core__stream_linenumber(T_sp tstream);
    T_sp cl__close(T_sp strm, T_sp abort);
    T_sp cl__open(T_sp filename,
             T_sp direction,
             T_sp element_type,
             T_sp if_exists, bool iesp,
             T_sp if_does_not_exist, bool idnesp,
             T_sp external_format,
             T_sp cstream);
    bool cl__streamp(T_sp strm);
    T_sp cl__interactive_stream_p(T_sp strm);
    T_sp cl__output_stream_p(T_sp strm);
    T_sp cl__input_stream_p(T_sp strm);
    T_sp cl__file_position(T_sp stream, T_sp position);
    T_sp cl__file_length(T_sp strm);
    T_sp cl__read_sequence(T_sp sequence, T_sp stream, T_sp start, T_sp oend);
    T_sp core__do_write_sequence(T_sp seq, T_sp stream, T_sp s, T_sp e);
    T_sp core__file_stream_fd(T_sp s);
    T_sp core__set_buffering_mode(T_sp stream, T_sp buffer_mode_symbol);
    T_sp cl__synonym_stream_symbol(T_sp strm);
    T_sp cl__make_synonym_stream(T_sp tsym);
    T_sp cl__make_concatenated_stream(List_sp ap);
    T_sp cl__make_echo_stream(T_sp strm1, T_sp strm2);
    T_sp cl__make_broadcast_stream(List_sp ap);
    T_sp cl__two_way_stream_output_stream(T_sp strm);
    T_sp cl__two_way_stream_input_stream(T_sp strm);
    T_sp cl__make_two_way_stream(T_sp istrm, T_sp ostrm);
    T_sp cl__make_string_input_stream(Str_sp strng, Fixnum_sp istart, T_sp iend);
    T_sp cl__get_output_stream_string(T_sp strm);
    T_sp cl__make_string_output_stream(Symbol_sp elementType);
    T_sp core__make_string_output_stream_from_string(T_sp s);
    T_mv core__process_lambda_list(List_sp lambdaList, T_sp context);
    T_mv core__process_single_dispatch_lambda_list(List_sp lambda_list);
    T_sp core__process_macro_lambda_list(List_sp lambda_list);
    Str_sp core__magic_name(const std::string& name);
    Str_sp core__lispify_name(Str_sp name);
    int core__environment_id(T_sp tenv);
    T_sp core__runtime_environment(T_sp tenv);
    T_sp core__environment_type_list(T_sp env);
    T_sp core__environment_list(T_sp env);
    T_sp core__environment_activation_frame(T_sp env);
    T_mv core__lexical_macro_function(T_sp name, T_sp env);
    T_mv core__lexical_function(T_sp name, T_sp env);
    T_sp core__environment_debug_values(T_sp frame);
    T_sp core__environment_debug_names(T_sp frame);
    int core__environment_length(T_sp frame);
    T_mv core__classify_return_from_symbol(T_sp env, Symbol_sp sym);
    T_sp cl__nreconc(List_sp list, T_sp tail);
    T_sp cl__revappend(List_sp list, T_sp tail);
    T_sp cl__nconc(List_sp lists);
    T_sp cl__last(T_sp list, int n);
    T_sp cl__listSTAR(T_sp tobjects);
    T_sp cl__list(T_sp objects);
    List_sp cl__nbutlast(List_sp l, Integer_sp in);
    List_sp cl__butlast(List_sp ll, Integer_sp in);
    T_sp cl__copy_list(T_sp arg);
    T_sp cl__nthcdr(int idx, T_sp arg);
    T_sp cl__nth(int idx, T_sp arg);
    T_sp cl__rassoc(T_sp item, List_sp a_list, T_sp test, T_sp test_not, T_sp key);
    T_sp core__function_source_code(Function_sp fn);
    T_mv cl__function_lambda_expression(Function_sp fn);
    void core__set_kind(Function_sp fn, Symbol_sp kind);
    gc::Nilable<SourcePosInfo_sp> core__function_source_pos_info(T_sp functionDesignator);
    T_mv core__function_lambda_list(T_sp obj);
    T_sp core__mkdir(T_sp directory, T_sp mode);
    Ratio_sp core__unix_get_local_time_zone();
    bool core__unix_daylight_saving_time(Integer_sp unix_time);
    T_sp cl__directory(T_sp mask, T_sp resolveSymlinks);
    T_sp core__copy_file(T_sp orig, T_sp dest);
    void core__chmod(T_sp file, T_sp mode);
    T_sp core__rmdir(T_sp directory);
    T_sp core__mkstemp(Str_sp thetemplate);
    Pathname_sp cl__user_homedir_pathname(T_sp host);
    T_sp cl__file_author(T_sp file);
    Number_sp cl__file_write_date(T_sp pathspec);
    T_sp cl__probe_file(T_sp filespec);
    T_sp cl__delete_file(T_sp file);
    T_mv cl__rename_file(T_sp oldn, T_sp newn, T_sp if_exists);
    Pathname_sp cl__truename(T_sp orig_pathname);
    Pathname_mv core__file_truename(T_sp pathname, T_sp filename, bool follow_links);
    T_sp core__readlink(Str_sp filename);
    Symbol_sp core__file_kind(T_sp filename, bool follow_links);
    Str_sp core__current_dir();
    T_sp ext__chdir(Pathname_sp dir);
    T_sp core__getppid();
    T_sp core__getpid();
    int core__waitpid(Fixnum_sp pid, Fixnum_sp options);
    T_sp core__fork();
    SourceManager_sp core__make_source_manager();
    void core__dump_source_manager(T_sp dumpAll);
    T_sp core__walk_to_find_source_pos_info(T_sp obj, T_sp defaultSpi);
    void core__walk_to_assign_source_pos_info(T_sp obj, SourcePosInfo_sp topInfo, T_sp stream);
    T_mv core__walk_to_find_source_info(T_sp obj);
    Fixnum_sp core__source_pos_info_column(SourcePosInfo_sp info);
    Fixnum_sp core__source_pos_info_lineno(SourcePosInfo_sp info);
    Integer_sp core__source_pos_info_filepos(SourcePosInfo_sp info);
    T_mv core__source_file_info(T_sp sourceFile, T_sp sourceDebugNamestring, size_t sourceDebugOffset, bool useLineno);
    T_sp cl__print(T_sp obj, T_sp output_stream_desig);
    T_sp cl__prin1(T_sp obj, T_sp output_stream_desig);
    T_sp cl__princ(T_sp obj, T_sp output_stream_desig);
    void cl__pprint(T_sp obj, T_sp stream);
    void core__print_unreadable_object_function(T_sp o, T_sp ostream, T_sp type, T_sp id, T_sp function);
    T_sp cl__write(T_sp x, T_sp strm, T_sp array, T_sp base,
              T_sp cas, T_sp circle, T_sp escape, T_sp gensym, T_sp length,
              T_sp level, T_sp lines, T_sp miser_width, T_sp pprint_dispatch,
              T_sp pretty, T_sp radix, T_sp readably, T_sp right_margin);
    T_sp core__write_object(T_sp obj, T_sp ostrm);
    T_sp cl__load(T_sp source, T_sp verbose, T_sp print, T_sp if_does_not_exist, T_sp external_format, T_sp search_list);
    T_sp core__load_source(T_sp source, bool verbose, bool print, T_sp externalFormat);
    Real_sp cl__imagpart(Number_sp x);
    Real_sp cl__realpart(Number_sp x);
    Complex_sp cl__complex(Real_sp r, Real_sp i);
    Real_mv cl__integer_decode_float(Float_sp x);
    Integer_sp cl__float_precision(Float_sp x);
    Integer_sp cl__float_digits(Float_sp x);
    Float_sp cl__float_sign(Float_sp x, Float_sp y, T_sp yp);
    Number_sp cl__scale_float(Number_sp x, Number_sp y);
    Number_mv cl__decode_float(Float_sp x);
    Real_sp cl__rem(Real_sp x, Real_sp y);
    Real_sp cl__mod(Real_sp x, Real_sp y);
    Number_mv cl__round(Real_sp x, T_sp y);
    Real_mv cl__truncate(Real_sp x, T_sp y);
    Real_mv cl__ceiling(Real_sp x, T_sp y);
    Real_mv cl__floor(Real_sp x, T_sp y);
    Number_sp cl__denominator(Number_sp x);
    Number_sp cl__numerator(Number_sp x);
    Float_sp cl__float(Real_sp x, T_sp y);
    StrWithFillPtr_sp core__integer_to_string(StrWithFillPtr_sp buffer, Integer_sp integer,
                                       Fixnum_sp base, bool radix, bool decimalp);
    StrWithFillPtr_sp core__bignum_to_string(StrWithFillPtr_sp buffer, const Bignum &bn, Fixnum_sp base);
    Integer_sp cl__lcm(List_sp nums);
    Integer_sp cl__gcd(List_sp nums);
    double core__num_op_atanh(double x);
    double core__num_op_acosh(double x);
    double core__num_op_asinh(double x);
    double core__num_op_acos(double x);
    double core__num_op_asin(double x);
    double randomNumberNormal01();
    double randomNumber01();
    void seedRandomNumberGeneratorsUsingTime();
    void seedRandomNumberGenerators(uint i);
    Integer_sp cl__get_universal_time();
    gc::Fixnum cl__integer_length(Integer_sp i);
    Number_sp core__log1p(Number_sp arg);
    Number_sp cl__log(Number_sp number, T_sp base);
    T_sp cl__atan(Number_sp x, T_sp y);
    Number_sp cl__expt(Number_sp x, Number_sp y);
    Number_sp cl__exp(Number_sp x);
    Number_sp cl__conjugate(Number_sp x);
    Number_sp cl__tanh(Number_sp x);
    Number_sp cl__cosh(Number_sp x);
    Number_sp cl__sinh(Number_sp x);
    Number_sp cl__tan(Number_sp x);
    Number_sp cl__cos(Number_sp x);
    Number_sp cl__sin(Number_sp x);
    Number_sp cl__sqrt(Number_sp x);
    DoubleFloat_mv core__nan();
    T_sp cl___EQ_(List_sp args);
    T_sp cl___NE_(List_sp args);
    T_mv cl___GE_(List_sp args);
    T_mv cl___LE_(List_sp args);
    T_mv cl___GT_(List_sp args);
    T_sp cl___LT_(List_sp args);
    bool two_arg__EQ_(Number_sp x, Number_sp y);
    bool two_arg__GE_(Number_sp x, Number_sp y);
    bool two_arg__GT_(Number_sp x, Number_sp y);
    bool two_arg__LE_(Number_sp x, Number_sp y);
    bool two_arg__LT_(Number_sp x, Number_sp y);
    T_sp cl___DIVIDE_(Number_sp num, List_sp numbers);
    T_mv cl___MINUS_(Number_sp num, List_sp numbers);
    T_mv cl___TIMES_(List_sp numbers);
    T_mv cl___PLUS_(List_sp numbers);
    Number_sp contagen_div(Number_sp na, Number_sp nb);
    Number_sp contagen_mul(Number_sp na, Number_sp nb);
    Number_sp contagen_sub(Number_sp na, Number_sp nb);
    Number_sp contagen_add(Number_sp na, Number_sp nb);
    T_mv cl__lognor(Integer_sp a, Integer_sp b);
    T_mv cl__lognand(Integer_sp a, Integer_sp b);
    T_mv cl__lognot(Integer_sp a);
    T_mv cl__logorc2(Integer_sp a, Integer_sp b);
    T_mv cl__logorc1(Integer_sp a, Integer_sp b);
    T_mv cl__logandc2(Integer_sp a, Integer_sp b);
    T_mv cl__logandc1(Integer_sp a, Integer_sp b);
    Integer_mv cl__logeqv(List_sp integers);
    Integer_sp cl__logxor(List_sp integers);
    Integer_sp cl__logior(List_sp integers);
    Integer_sp cl__logand(List_sp integers);
    Real_sp cl__max(Real_sp max, List_sp nums);
    Real_sp cl__min(Real_sp min, List_sp nums);
    Fixnum_sp core__fixnum_number_of_bits();
    Integer_sp core__convert_overflow_result_to_bignum(Fixnum_sp z);
    bool cl__zerop(T_sp num);
    void core__hash_table_force_rehash(HashTable_sp ht);
    T_mv cl__gethash(T_sp key, T_sp hashTable, T_sp default_value);
    bool cl__remhash(T_sp key, HashTable_sp ht);
    int core__hash_equalp(List_sp args);
    int core__hash_equal(List_sp args);
    int core__hash_eql(List_sp args);
    bool core__hash_table_entry_deleted_p(T_sp cons);
    T_mv cl__clrhash(HashTable_sp hash_table);
    T_mv cl__maphash(T_sp function_desig, T_sp thash_table);
    Symbol_sp core__hash_table_weakness(T_sp ht);
    T_sp cl__make_hash_table(T_sp test, Fixnum_sp size, Number_sp rehash_size, Real_sp orehash_threshold, Symbol_sp weakness, T_sp debug);
    void core__debug_hash_table(bool don);
    SmallMultimap_sp core__make_small_multimap();
    SmallMap_sp core__make_small_map();
    void core__print_current_ihs_frame_environment();
    void core__frame_pointers();
    void core__clib_backtrace(int depth);
    void core__low_level_backtrace();
    T_sp cl__random(T_sp olimit, RandomState_sp random_state);
    T_sp core__ihs_backtrace(T_sp outputDesignator, T_sp msg);
    void core__dynamic_binding_stack_dump(std::ostream &out);
    void core__exception_stack_dump();
    Vector_sp core__exception_stack();
    T_sp core__bds_val(int idx);
    Symbol_sp core__bds_var(int idx);
    int core__bds_top();
    int core__set_ihs_current_frame(int icf);
    int core__ihs_current_frame();
    int core__ihs_bds(int idx);
    T_sp core__ihs_env(int idx);
    T_sp core__ihs_arguments(int idx);
    T_sp core__ihs_fun(int idx);
    int core__ihs_next(int idx);
    int core__ihs_prev(int idx);
    int core__ihs_top();
    void core__ihs_backtrace_no_args();
    InvocationHistoryFrameIterator_sp core__get_invocation_history_frame_next(int idx);
    InvocationHistoryFrameIterator_sp core__get_invocation_history_frame_prev(int idx);
    InvocationHistoryFrameIterator_sp core__get_invocation_history_frame(int idx);
    InvocationHistoryFrameIterator_sp core__get_invocation_history_frame_top();
    InvocationHistoryFrameIterator_sp core__get_invocation_history_frame_search(T_sp idx, Symbol_sp direction);
    Integer_sp cl__sxhash(T_sp obj);
    T_sp cl__type_of(T_sp x);
    Symbol_mv core__type_to_symbol(T_sp x);
    Symbol_mv cl__gensym(T_sp x);
    T_mv core__sequence_start_end(T_sp func, T_sp sequence, Fixnum_sp start, T_sp end);
    List_sp cl__append(List_sp lists);
    T_mv cl__mapcan(T_sp op, List_sp lists);
    T_mv cl__mapcon(T_sp op, List_sp lists);
    T_mv core__mapappend(Function_sp fun, List_sp cargs);
    T_sp cl__mapl(T_sp top, List_sp lists);
    T_sp cl__maplist(T_sp func_desig, List_sp lists);
    T_sp cl__mapc(T_sp top, List_sp lists);
    T_sp cl__mapcar(T_sp func_desig, List_sp lists);
    T_sp cl__notevery(T_sp predicate, List_sp sequences);
    T_sp cl__notany(T_sp predicate, List_sp sequences);
    T_sp cl__some(T_sp predicate, List_sp sequences);
    T_sp cl__every(T_sp predicate, List_sp sequences);
    T_sp cl__read_preserving_whitespace(T_sp input_stream_designator, T_sp eof_error_p, T_sp eof_value, T_sp recursive_p);
    T_sp cl__read(T_sp input_stream_designator, T_sp eof_error_p, T_sp eof_value, T_sp recursive_p);
    T_mv cl__read_delimited_list(Character_sp chr, T_sp input_stream_designator, T_sp recursive_p);
    T_mv cl__fmakunbound(T_sp functionName);
    bool cl__fboundp(T_sp functionName);
    T_sp cl__fdefinition(T_sp functionName);
    T_sp core__STARfset(T_sp functionName, Function_sp functionObject, T_sp macro);
    Class_sp cl__class_of(T_sp obj);
    T_mv cl__null(T_sp obj);
    T_mv core__macroexpand_default(Function_sp macro_function, T_sp form, T_sp macro_env);
    T_mv cl__identity(T_sp arg);
    bool cl__constantp(T_sp obj, T_sp env);
    void core__gdb_inspect(Str_sp msg, T_sp o);
    void core__trap_execution(T_sp msg);
    Integer_sp core__cxx_lambda_list_handler_create_bindings_calls();
    void core__gdb(T_sp msg);
    void core__break(T_sp fmt, List_sp args);
    Integer_sp cl__ash(Integer_sp integer, Integer_sp count);
    T_mv cl__special_operator_p(T_sp sym);
    T_sp cl__macro_function(Symbol_sp symbol, T_sp env);
    T_mv core__separate_pair_list(List_sp listOfPairs);
    void core__test_memory_error();
    T_mv core__valid_function_name_p(T_sp arg);
    Symbol_mv core__function_block_name(T_sp functionName);
    T_mv cl__values_list(List_sp list);
    T_mv core__values_testing(List_sp args);
    T_mv cl__values(List_sp args);
    T_mv core__smart_pointer_details();
    T_sp core__unbound();
    T_mv core__substitute();
    bool core__is_true(T_sp arg);
    int core__pointer(T_sp obj);
    void core__setenv(Str_sp name, Str_sp arg, bool overwrite);
    void core__describe_cxx_object(T_sp obj, T_sp stream);
    T_sp ext__llvm_version();
    T_sp ext__getenv(Str_sp arg);
    T_mv ext__vfork_execvp(List_sp call_and_arguments);
    T_mv ext__system(Str_sp cmd);
    void ext__setenv(Str_sp name, Str_sp value);
    void core__dump_tagged_fixnum(T_sp val);
    gctools::Fixnum core__from_tagged_fixnum(T_sp val);
    T_sp core__to_tagged_fixnum(int val);
    Vector_sp core__all_registered_class_names();
    int core__incomplete_next_higher_power_of_2(Fixnum_sp fn);
    void core__dump_address_of(T_sp arg);
    T_sp cl__set(Symbol_sp sym, T_sp val);
    Str_sp core__argv(int idx);
    int core__argc();
    T_sp cl__machine_instance();
    T_sp cl__machine_version();
    T_sp cl__machine_type();
    T_sp cl__software_version();
    T_sp cl__software_type();
    T_sp core__create_tagged_immediate_value_or_nil(T_sp object);
    T_sp core__lisp_implementation_id();
    T_mv core__single_dispatch_method_cache_status();
    T_mv core__slot_cache_status();
    T_mv core__method_cache_status();
    void core__single_dispatch_method_cache_resize(Fixnum_sp pow);
    void core__slot_cache_resize(Fixnum_sp pow);
    void core__method_cache_resize(Fixnum_sp pow);
    T_sp cl__lisp_implementation_version();
    T_sp cl__lisp_implementation_type();
    void cl__sleep(T_sp oseconds);
    void core__clear_generic_function_dispatch_cache();
    Instance_sp core__copy_instance(Instance_sp obj);
    T_sp core__instance_class_set(T_sp obj, Class_sp mc);
    T_sp clos__setFuncallableInstanceFunction(T_sp obj, T_sp func);
    void core__clear_gfun_hash(T_sp what);
    T_sp core__maybe_expand_generic_function_arguments(T_sp args);
    T_mv core__progv_function(List_sp symbols, List_sp values, Function_sp func);
    void core__throw_function(T_sp tag, T_sp result_form);
    T_mv core__catch_function(T_sp tag, Function_sp thunk);
    T_mv core__multiple_value_prog1_function(Function_sp func1, Function_sp func2);
    T_mv core__multiple_value_funcall(T_sp funcDesignator, List_sp functions);
    T_mv core__funwind_protect(T_sp protected_fn, T_sp cleanup_fn);
    T_mv core__call_with_variable_bound(Symbol_sp sym, T_sp val, T_sp thunk);
    T_mv compiler__implicit_compile_hook_default(T_sp form, T_sp env);
    T_mv core__dladdr(Integer_sp addr);
    void core__call_dl_main_function(Pointer_sp addr);
    T_sp core__dlsym(T_sp ohandle, Str_sp name);
    T_mv core__dlopen(T_sp pathDesig);
    T_mv core__load_bundle(T_sp pathDesig, T_sp verbose, T_sp print, T_sp external_format);
    T_sp core__startup_image_pathname();
    T_mv core__mangle_name(Symbol_sp sym, bool is_function);
    Integer_sp core__cxx_fibn(Fixnum_sp reps, Fixnum_sp num);
    Fixnum_sp core__test_tagged_cast(Fixnum_sp pow2);
    void core__help_booting();
    T_mv core__eval_with_env_default(T_sp form, T_sp env);
    void core__evaluate_verbosity(Fixnum_sp level);
    int core__evaluate_depth();
    T_mv core__classify_let_variables_and_declares(List_sp variables, List_sp declaredSpecials);
    T_sp core__lookup_symbol_macro(Symbol_sp sym, T_sp env);
    T_sp core__extract_lambda_name(List_sp lambdaExpression, T_sp defaultValue);
    T_sp core__extract_lambda_name_from_declares(List_sp declares, T_sp defaultValue);
    T_mv core__process_declarations(List_sp inputBody, T_sp expectDocString);
    Function_sp core__coerce_to_function(T_sp arg);
    T_mv cl__funcall(T_sp function_desig, VaList_sp args);
    T_mv cl__eval(T_sp form);
    T_mv cl__apply(T_sp head, VaList_sp args);
    T_mv core__compile_form_and_eval_with_env(T_sp form, T_sp env, T_sp stepping, T_sp compiler_env_p, T_sp execute);
    void core__weak_rehash(WeakKeyHashTable_sp ht, T_sp sz);
    void core__weak_splat(WeakKeyHashTable_sp ht, Fixnum_sp idx);
    void core__weak_clrhash(WeakKeyHashTable_sp ht);
    void core__weak_remhash(WeakKeyHashTable_sp ht, T_sp key);
    void core__weak_setf_gethash(T_sp key, WeakKeyHashTable_sp ht, T_sp val);
    T_mv core__weak_gethash(T_sp tkey, WeakKeyHashTable_sp ht, T_sp defaultValue);
    WeakKeyHashTable_sp core__make_weak_key_hash_table(Fixnum_sp size);
    VectorDisplaced_sp core__make_vector_displaced(T_sp dim, T_sp elementType, T_sp displacedTo, size_t displacedOffset);
    int core__test_add(int x, int y);
    inline bool cl__listp(T_sp obj);
    inline bool cl__characterp(T_sp obj);
    inline bool cl__consp(T_sp obj);
    bool cl__atom(T_sp obj);
    inline bool cl__endp(T_sp arg);
    inline bool core__fixnump(T_sp arg);
    inline bool cl__symbolp(T_sp obj);
    inline Number_sp clasp_negate(Number_sp num);
    inline Number_sp clasp_one_minus(Number_sp num);
    inline Number_sp clasp_one_plus(Number_sp num);
    inline Number_sp clasp_signum(Number_sp num);
    inline Number_sp clasp_abs(Number_sp num);
    inline bool clasp_oddp(Integer_sp num);
    inline bool clasp_evenp(Integer_sp num);
    inline bool clasp_minusp(Real_sp num);
    inline bool clasp_plusp(Real_sp num);
    inline T_sp oTenth(T_sp o);
    inline T_sp oNinth(T_sp o);
    inline T_sp oEighth(T_sp o);
    inline T_sp oSeventh(T_sp o);
    inline T_sp oSixth(T_sp o);
    inline T_sp oFifth(T_sp o);
    inline T_sp oFourth(T_sp o);
    inline T_sp oThird(T_sp o);
    inline T_sp oSecond(T_sp o);
    inline T_sp oFirst(T_sp o);
    inline T_sp oCddddr(T_sp o);
    inline T_sp oCddadr(T_sp o);
    inline T_sp oCdaddr(T_sp o);
    inline T_sp oCdaadr(T_sp o);
    inline T_sp oCadddr(T_sp o);
    inline T_sp oCadadr(T_sp o);
    inline T_sp oCaaddr(T_sp o);
    inline T_sp oCaaadr(T_sp o);
    inline T_sp oCdddar(T_sp o);
    inline T_sp oCddaar(T_sp o);
    inline T_sp oCdadar(T_sp o);
    inline T_sp oCdaaar(T_sp o);
    inline T_sp oCaddar(T_sp o);
    inline T_sp oCadaar(T_sp o);
    inline T_sp oCaadar(T_sp o);
    inline T_sp oCaaaar(T_sp o);
    inline T_sp oCdddr(T_sp o);
    inline T_sp oCddar(T_sp o);
    inline T_sp oCdadr(T_sp o);
    inline T_sp oCdaar(T_sp o);
    inline T_sp oCaddr(T_sp o);
    inline T_sp oCadar(T_sp o);
    inline T_sp oCaadr(T_sp o);
    inline T_sp oCaaar(T_sp o);
    inline T_sp oCddr(T_sp o);
    inline T_sp oCdar(T_sp o);
    inline T_sp oCadr(T_sp o);
    inline T_sp oCaar(T_sp o);
    inline T_sp cl__rest(List_sp obj);
    inline T_sp oCdr(List_sp obj);
    inline T_sp oCar(List_sp obj);
    bool cl__equalp(T_sp x, T_sp y);
    bool cl__equal(T_sp x, T_sp y);
    bool cl__eql(T_sp x, T_sp y);
    bool cl__eq(T_sp x, T_sp y);
};
namespace ext {
    core::T_sp ext__maybeQuote(core::T_sp form);
};
namespace units {
};
namespace adapt {
    IndexedObjectBag_sp chem__makeIndexedObjectBag();
};
namespace geom {
};
namespace chem {
    core::T_sp chem__oligomer(Oligomer_O::NameType::smart_ptr oligomerName, core::List_sp parts);
    core::T_sp chem__calculate_point();
    core::T_sp chem__link(core::Symbol_sp monomer1Id, core::Symbol_sp coupling, OligomerPart_Monomer_sp monomer2 );
    core::T_sp chem__monomer(core::Symbol_sp monomerId, core::Symbol_sp groupName, core::List_sp monomerAliases, core::Str_sp comment);
    core::T_sp chem__atom_pos(core::List_sp args);
    core::T_sp chem__find_residue(core::List_sp args);
    core::T_sp chem__save_mol2(Matter_sp matter, core::T_sp destDesig);
    core::T_sp chem__load_mol2(core::T_sp fileName);
    geom::CoordinateArray_sp chem__make_coordinate_array_from_atom_list(core::List_sp atoms);
    core::T_sp chem__save_archive_with_auto_set_cando_database();
    core::T_sp chem__create_virtual_atom();
    core::T_sp chem__bundle_database_path(core::T_sp pathDesig);
    core::T_sp chem__standard_database(core::T_sp pathDesig, core::T_sp loadSeed, core::Fixnum_sp overbosity);
    core::T_sp chem__database(core::T_sp fileNameDesig, core::Fixnum_sp overbosity);
    ReadAmberParameters_sp chem__make_read_amber_parameters();
    core::T_sp chem__set_monomer_pack(core::Symbol_sp packName, core::List_sp parts, core::Cons_sp atomNames );
    core::T_sp chem__extend_aliases(core::Symbol_sp packName, core::List_sp parts, core::Cons_sp atomAliases);
    core::T_sp chem__define_monomer_pack(core::Symbol_sp packName, core::List_sp parts, core::Cons_sp atomAliases );
    core::List_sp chem__parseChimeraAtomSpecs(const string& spec);
    core::T_sp chem__alias(core::Symbol_sp monAlias, core::Symbol_sp atomAlias);
    core::T_sp chem__create_for_matter(Matter_sp matter);
    core::T_sp chem__load_pdb(core::Str_sp fileName);
    core::T_sp chem__set_oligomer(Oligomer_O::NameType::smart_ptr oligomerName, core::List_sp parts);
    core::T_sp chem__oligomer_sequence(Oligomer_sp olig);
    adapt::ObjectSet_sp chem__atomsInMatterAsObjectSet( Matter_sp matter);
    adapt::ObjectSet_sp chem__atomsWithinSphereAsObjectSet( Matter_sp matter, Vector3 center, double radius );
    void chem__connectAtomsInMatterInCovalentContact(Matter_sp matter);
    bool chem__is_in_plug_name(core::Symbol_sp plugName);
    core::T_sp chem__map_atoms(core::Symbol_sp result_type, core::T_sp funcDesig, Matter_sp m);
    core::T_sp chem__map_residues(core::Symbol_sp result_type, core::T_sp funcDesig, Matter_sp m);
    core::T_sp chem__map_molecules(core::Symbol_sp result_type, core::T_sp funcDesig, Matter_sp m );
    PdbMonomerDatabase_sp chem__readPdbMonomerConnectivityDatabase(const string& fileName);
    Element elementFromAtomNameString(const string& name);
    Element elementFromAtomNameStringCaseInsensitive(const string& name);
    double vdwRadiusForElement(Element element);
    FFTypesDb_sp chem__make_fftypes_db();
};
#endif // EXPOSE_FUNCTION_SIGNATURES
#ifdef EXPOSE_FUNCTION_BINDINGS_HELPERS
NOINLINE void expose_function_1_helper() {
  expose_function(core::magic_name("ast_tooling__testDerivable"),true,&asttooling::ast_tooling__testDerivable,"");
}
NOINLINE void expose_function_2_helper() {
  expose_function(core::magic_name("ast_tooling__deduplicate"),true,&asttooling::ast_tooling__deduplicate,"");
}
NOINLINE void expose_function_3_helper() {
  expose_function(core::magic_name("ast_tooling__clangVersionString"),true,&asttooling::ast_tooling__clangVersionString,"");
}
NOINLINE void expose_function_4_helper() {
  expose_function(core::magic_name("ast_tooling__makeAstVisitor"),true,&asttooling::ast_tooling__makeAstVisitor,"");
}
NOINLINE void expose_function_5_helper() {
  expose_function(core::magic_name("serve_event_internal__ll_serveEventWithTimeout"),true,&serveEvent::serve_event_internal__ll_serveEventWithTimeout,"");
}
NOINLINE void expose_function_6_helper() {
  expose_function(core::magic_name("serve_event_internal__ll_serveEventNoTimeout"),true,&serveEvent::serve_event_internal__ll_serveEventNoTimeout,"");
}
NOINLINE void expose_function_7_helper() {
  expose_function(core::magic_name("serve_event_internal__ll_fdset_size"),true,&serveEvent::serve_event_internal__ll_fdset_size,"");
}
NOINLINE void expose_function_8_helper() {
  expose_function(core::magic_name("serve_event_internal__ll_fd_isset"),true,&serveEvent::serve_event_internal__ll_fd_isset,"");
}
NOINLINE void expose_function_9_helper() {
  expose_function(core::magic_name("serve_event_internal__ll_fd_set"),true,&serveEvent::serve_event_internal__ll_fd_set,"");
}
NOINLINE void expose_function_10_helper() {
  expose_function(core::magic_name("serve_event_internal__ll_fd_zero"),true,&serveEvent::serve_event_internal__ll_fd_zero,"");
}
NOINLINE void expose_function_11_helper() {
  expose_function(core::magic_name("sockets_internal__ll_setSockoptLinger"),true,&sockets::sockets_internal__ll_setSockoptLinger,"(fd level constant value)");
}
NOINLINE void expose_function_12_helper() {
  expose_function(core::magic_name("sockets_internal__ll_setSockoptTimeval"),true,&sockets::sockets_internal__ll_setSockoptTimeval,"(fd level constant value)");
}
NOINLINE void expose_function_13_helper() {
  expose_function(core::magic_name("sockets_internal__ll_setSockoptBool"),true,&sockets::sockets_internal__ll_setSockoptBool,"(fd level constant value)");
}
NOINLINE void expose_function_14_helper() {
  expose_function(core::magic_name("sockets_internal__ll_setSockoptInt"),true,&sockets::sockets_internal__ll_setSockoptInt,"(fd level constant value)");
}
NOINLINE void expose_function_15_helper() {
  expose_function(core::magic_name("sockets_internal__ll_getSockoptLinger"),true,&sockets::sockets_internal__ll_getSockoptLinger,"(fd level constant)");
}
NOINLINE void expose_function_16_helper() {
  expose_function(core::magic_name("sockets_internal__ll_getSockoptTimeval"),true,&sockets::sockets_internal__ll_getSockoptTimeval,"(fd level constant)");
}
NOINLINE void expose_function_17_helper() {
  expose_function(core::magic_name("sockets_internal__ll_getSockoptBool"),true,&sockets::sockets_internal__ll_getSockoptBool,"(fd level constant)");
}
NOINLINE void expose_function_18_helper() {
  expose_function(core::magic_name("sockets_internal__ll_getSockoptInt"),true,&sockets::sockets_internal__ll_getSockoptInt,"(fd level constant)");
}
NOINLINE void expose_function_19_helper() {
  expose_function(core::magic_name("sockets_internal__ll_strerror_errno"),true,&sockets::sockets_internal__ll_strerror_errno,"");
}
NOINLINE void expose_function_20_helper() {
  expose_function(core::magic_name("sockets_internal__ll_strerror"),true,&sockets::sockets_internal__ll_strerror,"(num)");
}
NOINLINE void expose_function_21_helper() {
  expose_function(core::magic_name("sockets_internal__ll_autoCloseTwoWayStream"),true,&sockets::sockets_internal__ll_autoCloseTwoWayStream,"(stream)");
}
NOINLINE void expose_function_22_helper() {
  expose_function(core::magic_name("sockets_internal__ll_makeStreamFromFd"),true,&sockets::sockets_internal__ll_makeStreamFromFd,"(name fd stream-mode element-type external-format)");
}
NOINLINE void expose_function_23_helper() {
  expose_function(core::magic_name("sockets_internal__ll_dup"),true,&sockets::sockets_internal__ll_dup,"(fd)");
}
NOINLINE void expose_function_24_helper() {
  expose_function(core::magic_name("sockets_internal__ll_setfNonBlockingMode"),true,&sockets::sockets_internal__ll_setfNonBlockingMode,"(fd nblock)");
}
NOINLINE void expose_function_25_helper() {
  expose_function(core::magic_name("sockets_internal__ll_nonBlockingMode"),true,&sockets::sockets_internal__ll_nonBlockingMode,"(fd)");
}
NOINLINE void expose_function_26_helper() {
  expose_function(core::magic_name("sockets_internal__socketPeername_localSocket"),true,&sockets::sockets_internal__socketPeername_localSocket,"(fd)");
}
NOINLINE void expose_function_27_helper() {
  expose_function(core::magic_name("sockets_internal__ll_socketConnect_localSocket"),true,&sockets::sockets_internal__ll_socketConnect_localSocket,"(fd family path)");
}
NOINLINE void expose_function_28_helper() {
  expose_function(core::magic_name("sockets_internal__ll_socketAccept_localSocket"),true,&sockets::sockets_internal__ll_socketAccept_localSocket,"(socket-file-descriptor)");
}
NOINLINE void expose_function_29_helper() {
  expose_function(core::magic_name("sockets_internal__ll_socketBind_localSocket"),true,&sockets::sockets_internal__ll_socketBind_localSocket,"(fd name family)");
}
NOINLINE void expose_function_30_helper() {
  expose_function(core::magic_name("sockets_internal__socketSendNoAddress"),true,&sockets::sockets_internal__socketSendNoAddress,"(fb buffer length oob eor dontroute dontwait nosignal confirm)");
}
NOINLINE void expose_function_31_helper() {
  expose_function(core::magic_name("sockets_internal__ll_socketSendAddress"),true,&sockets::sockets_internal__ll_socketSendAddress,"(arg)");
}
NOINLINE void expose_function_32_helper() {
  expose_function(core::magic_name("sockets_internal__ll_socketName"),true,&sockets::sockets_internal__ll_socketName,"(fd vector)");
}
NOINLINE void expose_function_33_helper() {
  expose_function(core::magic_name("sockets_internal__ll_socketPeername_inetSocket"),true,&sockets::sockets_internal__ll_socketPeername_inetSocket,"(fd vector)");
}
NOINLINE void expose_function_34_helper() {
  expose_function(core::magic_name("sockets_internal__ll_socketConnect_inetSocket"),true,&sockets::sockets_internal__ll_socketConnect_inetSocket,"(port ip0 ip1 ip2 ip3 socket-file-descriptor)");
}
NOINLINE void expose_function_35_helper() {
  expose_function(core::magic_name("sockets_internal__ll_socketAccept_inetSocket"),true,&sockets::sockets_internal__ll_socketAccept_inetSocket,"(sfd)");
}
NOINLINE void expose_function_36_helper() {
  expose_function(core::magic_name("sockets_internal__ll_socketBind_inetSocket"),true,&sockets::sockets_internal__ll_socketBind_inetSocket,"(port ip0 ip1 ip2 ip3 socketFileDescriptor)");
}
NOINLINE void expose_function_37_helper() {
  expose_function(core::magic_name("sockets_internal__ll_getProtocolByName"),true,&sockets::sockets_internal__ll_getProtocolByName,"(name)");
}
NOINLINE void expose_function_38_helper() {
  expose_function(core::magic_name("sockets_internal__ll_socketReceive"),true,&sockets::sockets_internal__ll_socketReceive,"(fd buffer length oob peek waitall)");
}
NOINLINE void expose_function_39_helper() {
  expose_function(core::magic_name("sockets_internal__ll_getHostByAddress"),true,&sockets::sockets_internal__ll_getHostByAddress,"(address host-ent stf_host_ent_name setf_host_ent_aliases setf_host_ent_address_type setf_host_ent_addresses)");
}
NOINLINE void expose_function_40_helper() {
  expose_function(core::magic_name("sockets_internal__ll_getHostByName"),true,&sockets::sockets_internal__ll_getHostByName,"(host-name host-ent setf_host_ent_name setf_host_ent_aliases setf_host_ent_address_type setf_host_ent_addresses)");
}
NOINLINE void expose_function_41_helper() {
  expose_function(core::magic_name("sockets_internal__ll_getNameServiceErrorMessage"),true,&sockets::sockets_internal__ll_getNameServiceErrorMessage,"(num)");
}
NOINLINE void expose_function_42_helper() {
  expose_function(core::magic_name("sockets_internal__ll_get_name_service_h_errno"),true,&sockets::sockets_internal__ll_get_name_service_h_errno,"");
}
NOINLINE void expose_function_43_helper() {
  expose_function(core::magic_name("sockets_internal__ll_socket_errno"),true,&sockets::sockets_internal__ll_socket_errno,"");
}
NOINLINE void expose_function_44_helper() {
  expose_function(core::magic_name("sockets_internal__funcName"),true,&sockets::sockets_internal__funcName,"(arg)");
}
NOINLINE void expose_function_45_helper() {
  expose_function(core::magic_name("sockets_internal__ff_close"),true,&sockets::sockets_internal__ff_close,"(filedes)");
}
NOINLINE void expose_function_46_helper() {
  expose_function(core::magic_name("sockets_internal__ff_listen"),true,&sockets::sockets_internal__ff_listen,"(sockfd backlog)");
}
NOINLINE void expose_function_47_helper() {
  expose_function(core::magic_name("sockets_internal__ff_socket"),true,&sockets::sockets_internal__ff_socket,"(domain type protocol)");
}
NOINLINE void expose_function_48_helper() {
  expose_function(core::magic_name("gctools__next_header_kind"),true,&gctools::gctools__next_header_kind,"");
}
NOINLINE void expose_function_49_helper() {
  expose_function(core::magic_name("gctools__debug_allocations"),true,&gctools::gctools__debug_allocations,"");
}
NOINLINE void expose_function_50_helper() {
  expose_function(core::magic_name("gctools__cleanup"),true,&gctools::gctools__cleanup,"");
}
NOINLINE void expose_function_51_helper() {
  expose_function(core::magic_name("gctools__garbage_collect"),true,&gctools::gctools__garbage_collect,"");
}
NOINLINE void expose_function_52_helper() {
  expose_function(core::magic_name("gctools__stack_depth"),true,&gctools::gctools__stack_depth,"");
}
NOINLINE void expose_function_53_helper() {
  expose_function(core::magic_name("cl__room"),true,&gctools::cl__room,"(&optional x (marker 0) msg)");
}
NOINLINE void expose_function_54_helper() {
  expose_function(core::magic_name("gctools__alloc_pattern_end"),true,&gctools::gctools__alloc_pattern_end,"");
}
NOINLINE void expose_function_55_helper() {
  expose_function(core::magic_name("gctools__alloc_pattern_begin"),true,&gctools::gctools__alloc_pattern_begin,"");
}
NOINLINE void expose_function_56_helper() {
  expose_function(core::magic_name("gctools__gc_marker"),true,&gctools::gctools__gc_marker,"(&optional marker)");
}
NOINLINE void expose_function_57_helper() {
  expose_function(core::magic_name("gctools__monitor_allocations"),true,&gctools::gctools__monitor_allocations,"(on &key (backtrace-start 0) (backtrace-count 0) (backtrace-depth 6))");
}
NOINLINE void expose_function_58_helper() {
  expose_function(core::magic_name("gctools__gc_info"),true,&gctools::gctools__gc_info,"(&optional x (marker 0))");
}
NOINLINE void expose_function_59_helper() {
  expose_function(core::magic_name("core__hardwired_kinds"),true,&gctools::core__hardwired_kinds,"");
}
NOINLINE void expose_function_60_helper() {
  expose_function(core::magic_name("core__inherits_from_instance"),true,&gctools::core__inherits_from_instance,"(obj)");
}
NOINLINE void expose_function_61_helper() {
  expose_function(core::magic_name("core__header_kind"),true,&gctools::core__header_kind,"");
}
NOINLINE void expose_function_62_helper() {
  expose_function(core::magic_name("gctools__bytes_allocated"),true,&gctools::gctools__bytes_allocated,"");
}
NOINLINE void expose_function_63_helper() {
  expose_function(core::magic_name("gctools__deallocate_unmanaged_instance"),true,&gctools::gctools__deallocate_unmanaged_instance,"");
}
NOINLINE void expose_function_64_helper() {
  expose_function(core::magic_name("gctools__bootstrap_kind_p"),true,&gctools::gctools__bootstrap_kind_p,"");
}
NOINLINE void expose_function_65_helper() {
  expose_function(core::magic_name("gctools__bootstrap_kind_symbols"),true,&gctools::gctools__bootstrap_kind_symbols,"");
}
NOINLINE void expose_function_66_helper() {
  expose_function(core::magic_name("gctools__max_bootstrap_kinds"),true,&gctools::gctools__max_bootstrap_kinds,"");
}
NOINLINE void expose_function_67_helper() {
  expose_function(core::magic_name("core__telemetry_count"),true,&telemetry::core__telemetry_count,"");
}
NOINLINE void expose_function_68_helper() {
  expose_function(core::magic_name("core__telemetry_dump"),true,&telemetry::core__telemetry_dump,"(&optional (begin 0) end)");
}
NOINLINE void expose_function_69_helper() {
  expose_function(core::magic_name("core__telemetry_labels"),true,&telemetry::core__telemetry_labels,"");
}
NOINLINE void expose_function_70_helper() {
  expose_function(core::magic_name("core__telemetry_follow"),true,&telemetry::core__telemetry_follow,"(address)");
}
NOINLINE void expose_function_71_helper() {
  expose_function(core::magic_name("core__telemetry_search_labels"),true,&telemetry::core__telemetry_search_labels,"(label &optional (begin 0) end)");
}
NOINLINE void expose_function_72_helper() {
  expose_function(core::magic_name("core__telemetry_search"),true,&telemetry::core__telemetry_search,"(addresses)");
}
NOINLINE void expose_function_73_helper() {
  expose_function(core::magic_name("core__telemetry_open"),true,&telemetry::core__telemetry_open,"(pathname)");
}
NOINLINE void expose_function_74_helper() {
  expose_function(core::magic_name("llvm_sys__viewCFG"),true,&llvmo::llvm_sys__viewCFG,"(fn &optional only)");
}
NOINLINE void expose_function_75_helper() {
  expose_function(core::magic_name("llvm_sys__disassembleSTAR"),true,&llvmo::llvm_sys__disassembleSTAR,"");
}
NOINLINE void expose_function_76_helper() {
  expose_function(core::magic_name("llvm_sys__getOrCreateExternalGlobal"),true,&llvmo::llvm_sys__getOrCreateExternalGlobal,"");
}
NOINLINE void expose_function_77_helper() {
  expose_function(core::magic_name("llvm_sys__throwIfMismatchedStructureSizes"),true,&llvmo::llvm_sys__throwIfMismatchedStructureSizes,"(&key tsp tmv ihf)");
}
NOINLINE void expose_function_78_helper() {
  expose_function(core::magic_name("llvm_sys__cxxDataStructuresInfo"),true,&llvmo::llvm_sys__cxxDataStructuresInfo,"");
}
NOINLINE void expose_function_79_helper() {
  expose_function(core::magic_name("llvm_sys__mangleSymbolName"),true,&llvmo::llvm_sys__mangleSymbolName,"");
}
NOINLINE void expose_function_80_helper() {
  expose_function(core::magic_name("LLVM-SYS:createBasicAliasAnalysisPass"),true,&llvm::createBasicAliasAnalysisPass,"");
}
NOINLINE void expose_function_81_helper() {
  expose_function(core::magic_name("LLVM-SYS:createTypeBasedAliasAnalysisPass"),true,&llvm::createTypeBasedAliasAnalysisPass,"");
}
NOINLINE void expose_function_82_helper() {
  expose_function(core::magic_name("LLVM-SYS:createLowerExpectIntrinsicPass"),true,&llvm::createLowerExpectIntrinsicPass,"");
}
NOINLINE void expose_function_83_helper() {
  expose_function(core::magic_name("LLVM-SYS:createEarlyCSEPass"),true,&llvm::createEarlyCSEPass,"");
}
NOINLINE void expose_function_84_helper() {
  expose_function(core::magic_name("LLVM-SYS:createCorrelatedValuePropagationPass"),true,&llvm::createCorrelatedValuePropagationPass,"");
}
NOINLINE void expose_function_85_helper() {
  expose_function(core::magic_name("LLVM-SYS:createVerifierPass"),true,&llvm::createVerifierPass,"");
}
NOINLINE void expose_function_86_helper() {
  expose_function(core::magic_name("LLVM-SYS:createConstantPropagationPass"),true,&llvm::createConstantPropagationPass,"");
}
NOINLINE void expose_function_87_helper() {
  expose_function(core::magic_name("LLVM-SYS:createTailCallEliminationPass"),true,&llvm::createTailCallEliminationPass,"");
}
NOINLINE void expose_function_88_helper() {
  expose_function(core::magic_name("LLVM-SYS:createScalarReplAggregatesPass"),true,&llvm::createScalarReplAggregatesPass,"");
}
NOINLINE void expose_function_89_helper() {
  expose_function(core::magic_name("LLVM-SYS:createSCCPPass"),true,&llvm::createSCCPPass,"");
}
NOINLINE void expose_function_90_helper() {
  expose_function(core::magic_name("LLVM-SYS:createReassociatePass"),true,&llvm::createReassociatePass,"");
}
NOINLINE void expose_function_91_helper() {
  expose_function(core::magic_name("LLVM-SYS:createPromoteMemoryToRegisterPass"),true,&llvm::createPromoteMemoryToRegisterPass,"");
}
NOINLINE void expose_function_92_helper() {
  expose_function(core::magic_name("LLVM-SYS:createMemCpyOptPass"),true,&llvm::createMemCpyOptPass,"");
}
NOINLINE void expose_function_93_helper() {
  expose_function(core::magic_name("LLVM-SYS:createLoopUnswitchPass"),true,&llvm::createLoopUnswitchPass,"");
}
NOINLINE void expose_function_94_helper() {
  expose_function(core::magic_name("LLVM-SYS:createLoopUnrollPass"),true,&llvm::createLoopUnrollPass,"");
}
NOINLINE void expose_function_95_helper() {
  expose_function(core::magic_name("LLVM-SYS:createLoopRotatePass"),true,&llvm::createLoopRotatePass,"");
}
NOINLINE void expose_function_96_helper() {
  expose_function(core::magic_name("LLVM-SYS:createLoopIdiomPass"),true,&llvm::createLoopIdiomPass,"");
}
NOINLINE void expose_function_97_helper() {
  expose_function(core::magic_name("LLVM-SYS:createLoopDeletionPass"),true,&llvm::createLoopDeletionPass,"");
}
NOINLINE void expose_function_98_helper() {
  expose_function(core::magic_name("LLVM-SYS:createLICMPass"),true,&llvm::createLICMPass,"");
}
NOINLINE void expose_function_99_helper() {
  expose_function(core::magic_name("LLVM-SYS:createJumpThreadingPass"),true,&llvm::createJumpThreadingPass,"");
}
NOINLINE void expose_function_100_helper() {
  expose_function(core::magic_name("LLVM-SYS:createInstructionCombiningPass"),true,&llvm::createInstructionCombiningPass,"");
}
NOINLINE void expose_function_101_helper() {
  expose_function(core::magic_name("LLVM-SYS:createIndVarSimplifyPass"),true,&llvm::createIndVarSimplifyPass,"");
}
NOINLINE void expose_function_102_helper() {
  expose_function(core::magic_name("LLVM-SYS:createGVNPass"),true,&llvm::createGVNPass,"");
}
NOINLINE void expose_function_103_helper() {
  expose_function(core::magic_name("LLVM-SYS:createDeadStoreEliminationPass"),true,&llvm::createDeadStoreEliminationPass,"");
}
NOINLINE void expose_function_104_helper() {
  expose_function(core::magic_name("LLVM-SYS:createCFGSimplificationPass"),true,&llvm::createCFGSimplificationPass,"");
}
NOINLINE void expose_function_105_helper() {
  expose_function(core::magic_name("LLVM-SYS:createAggressiveDCEPass"),true,&llvm::createAggressiveDCEPass,"");
}
NOINLINE void expose_function_106_helper() {
  expose_function(core::magic_name("LLVM-SYS:InitializeNativeTarget"),true,&llvm::InitializeNativeTarget,"");
}
NOINLINE void expose_function_107_helper() {
  expose_function(core::magic_name("LLVM-SYS:createMemDepPrinter"),true,&llvm::createMemDepPrinter,"");
}
NOINLINE void expose_function_108_helper() {
  expose_function(core::magic_name("LLVM-SYS:createModuleDebugInfoPrinterPass"),true,&llvm::createModuleDebugInfoPrinterPass,"");
}
NOINLINE void expose_function_109_helper() {
  expose_function(core::magic_name("LLVM-SYS:createRegionInfoPass"),true,&llvm::createRegionInfoPass,"");
}
NOINLINE void expose_function_110_helper() {
  expose_function(core::magic_name("LLVM-SYS:createInstCountPass"),true,&llvm::createInstCountPass,"");
}
NOINLINE void expose_function_111_helper() {
  expose_function(core::magic_name("LLVM-SYS:createLazyValueInfoPass"),true,&llvm::createLazyValueInfoPass,"");
}
NOINLINE void expose_function_112_helper() {
  expose_function(core::magic_name("LLVM-SYS:createScalarEvolutionAliasAnalysisPass"),true,&llvm::createScalarEvolutionAliasAnalysisPass,"");
}
NOINLINE void expose_function_113_helper() {
  expose_function(core::magic_name("LLVM-SYS:createAAEvalPass"),true,&llvm::createAAEvalPass,"");
}
NOINLINE void expose_function_114_helper() {
  expose_function(core::magic_name("LLVM-SYS:createAlwaysInlinerPass"),true,(llvm::Pass * (*)()) & llvm::createAlwaysInlinerPass,"");
}
NOINLINE void expose_function_115_helper() {
  expose_function(core::magic_name("LLVM-SYS:createFunctionInliningPass"),true,(llvm::Pass * (*)(unsigned, unsigned)) & llvm::createFunctionInliningPass,"");
}
NOINLINE void expose_function_116_helper() {
  expose_function(core::magic_name("LLVM-SYS:createAliasAnalysisCounterPass"),true,&llvm::createAliasAnalysisCounterPass,"");
}
NOINLINE void expose_function_117_helper() {
  expose_function(core::magic_name("LLVM-SYS:TargetRegistryLookupTarget"),true,&llvmo::TargetRegistryLookupTarget,"");
}
NOINLINE void expose_function_118_helper() {
  expose_function(core::magic_name("LLVM-SYS:finalizeClosure"),true,&llvmo::finalizeClosure,"");
}
NOINLINE void expose_function_119_helper() {
  expose_function(core::magic_name("LLVM-SYS:finalizeEngineAndRegisterWithGcAndGetCompiledFunction"),true,&llvmo::finalizeEngineAndRegisterWithGcAndGetCompiledFunction,"");
}
NOINLINE void expose_function_120_helper() {
  expose_function(core::magic_name("LLVM-SYS:pointer-type-get"),true,&llvmo::PointerType_O::get,"(element-type &optional (address-space 0))");
}
NOINLINE void expose_function_121_helper() {
  expose_function(core::magic_name("LLVM-SYS:array-type-get"),true,&llvmo::ArrayType_O::get,"(element-type num-elements)");
}
NOINLINE void expose_function_122_helper() {
  expose_function(core::magic_name("LLVM-SYS:struct-type-get"),true,&llvmo::StructType_O::get,"");
}
NOINLINE void expose_function_123_helper() {
  expose_function(core::magic_name("LLVM-SYS:struct-type-create"),true,&llvmo::StructType_O::make,"(context &key elements name is-packed)");
}
NOINLINE void expose_function_124_helper() {
  expose_function(core::magic_name("LLVM-SYS:function-type-get"),true,&llvmo::FunctionType_O::get,"(result &optional params is-var-arg)");
}
NOINLINE void expose_function_125_helper() {
  expose_function(core::magic_name("LLVM-SYS:type-get-int64-ptr-ty"),true,&llvm::Type::getInt64PtrTy,"");
}
NOINLINE void expose_function_126_helper() {
  expose_function(core::magic_name("LLVM-SYS:type-get-int32-ptr-ty"),true,&llvm::Type::getInt32PtrTy,"");
}
NOINLINE void expose_function_127_helper() {
  expose_function(core::magic_name("LLVM-SYS:type-get-int8-ptr-ty"),true,&llvm::Type::getInt8PtrTy,"");
}
NOINLINE void expose_function_128_helper() {
  expose_function(core::magic_name("LLVM-SYS:type-get-int64-ty"),true,&llvm::Type::getInt64Ty,"");
}
NOINLINE void expose_function_129_helper() {
  expose_function(core::magic_name("LLVM-SYS:type-get-int32-ty"),true,&llvm::Type::getInt32Ty,"");
}
NOINLINE void expose_function_130_helper() {
  expose_function(core::magic_name("LLVM-SYS:type-get-int8-ty"),true,&llvm::Type::getInt8Ty,"");
}
NOINLINE void expose_function_131_helper() {
  expose_function(core::magic_name("LLVM-SYS:type-get-int1-ty"),true,&llvm::Type::getInt1Ty,"");
}
NOINLINE void expose_function_132_helper() {
  expose_function(core::magic_name("LLVM-SYS:type-get-void-ty"),true,&llvm::Type::getVoidTy,"");
}
NOINLINE void expose_function_133_helper() {
  expose_function(core::magic_name("LLVM-SYS:type-get-double-ty"),true,&llvm::Type::getDoubleTy,"");
}
NOINLINE void expose_function_134_helper() {
  expose_function(core::magic_name("LLVM-SYS:type-get-float-ty"),true,&llvm::Type::getFloatTy,"");
}
NOINLINE void expose_function_135_helper() {
  expose_function(core::magic_name("type-get-float-ty","LLVM-SYS"),true,&llvm::Type::getFloatTy,"");
}
NOINLINE void expose_function_136_helper() {
  expose_function(core::magic_name("LLVM-SYS:basic-block-create"),true,&llvm::BasicBlock::Create,"(context &optional (name \"\") parent basic_block)");
}
NOINLINE void expose_function_137_helper() {
  expose_function(core::magic_name("llvm_sys__FunctionCreate"),true,&llvmo::llvm_sys__FunctionCreate,"");
}
NOINLINE void expose_function_138_helper() {
  expose_function(core::magic_name("LLVM-SYS:value-as-metadata-get"),true,&llvmo::ValueAsMetadata_O::get,"");
}
NOINLINE void expose_function_139_helper() {
  expose_function(core::magic_name("LLVM-SYS:mdstring-get"),true,&llvmo::MDString_O::get,"");
}
NOINLINE void expose_function_140_helper() {
  expose_function(core::magic_name("LLVM-SYS:mdnode-get"),true,&llvmo::MDNode_O::get,"");
}
NOINLINE void expose_function_141_helper() {
  expose_function(core::magic_name("LLVM-SYS:make-irbuilder"),true,&llvmo::IRBuilder_O::make,"");
}
NOINLINE void expose_function_142_helper() {
  expose_function(core::magic_name("LLVM-SYS:makeAPInt64"),true,&llvmo::APInt_O::makeAPInt64,"");
}
NOINLINE void expose_function_143_helper() {
  expose_function(core::magic_name("LLVM-SYS:makeAPInt32"),true,&llvmo::APInt_O::makeAPInt32,"");
}
NOINLINE void expose_function_144_helper() {
  expose_function(core::magic_name("LLVM-SYS:makeAPIntWidth"),true,&llvmo::APInt_O::makeAPIntWidth,"");
}
NOINLINE void expose_function_145_helper() {
  expose_function(core::magic_name("LLVM-SYS:makeAPInt1"),true,&llvmo::APInt_O::makeAPInt1,"");
}
NOINLINE void expose_function_146_helper() {
  expose_function(core::magic_name("LLVM-SYS:makeAPInt"),true,&llvmo::APInt_O::makeAPInt,"");
}
NOINLINE void expose_function_147_helper() {
  expose_function(core::magic_name("LLVM-SYS:makeAPFloatDouble"),true,&llvmo::APFloat_O::makeAPFloatDouble,"(value)");
}
NOINLINE void expose_function_148_helper() {
  expose_function(core::magic_name("LLVM-SYS:make-apfloat-float"),true,&llvmo::APFloat_O::makeAPFloatFloat,"(value)");
}
NOINLINE void expose_function_149_helper() {
  expose_function(core::magic_name("LLVM-SYS:constant-pointer-null-get"),true,&llvm::ConstantPointerNull::get,"");
}
NOINLINE void expose_function_150_helper() {
  expose_function(core::magic_name("LLVM-SYS:UNDEF_VALUE-GET"),true,&llvm::UndefValue::get,"");
}
NOINLINE void expose_function_151_helper() {
  expose_function(core::magic_name("LLVM-SYS:CONSTANT-STRUCT-GET"),true,(llvm::Constant *(*)(llvm::StructType *T, llvm::ArrayRef<llvm::Constant *>)) &llvm::ConstantStruct::get,"");
}
NOINLINE void expose_function_152_helper() {
  expose_function(core::magic_name("LLVM-SYS:constant-int-get"),true,(llvm::ConstantInt *(*)(llvm::LLVMContext &, const llvm::APInt &)) &llvm::ConstantInt::get,"");
}
NOINLINE void expose_function_153_helper() {
  expose_function(core::magic_name("LLVM-SYS:constantFpGetTypeStringref"),true,(llvm::Constant *(*)(llvm::Type *, llvm::StringRef))&llvm::ConstantFP::get,"");
}
NOINLINE void expose_function_154_helper() {
  expose_function(core::magic_name("LLVM-SYS:constantFpGetTypeDouble"),true,(llvm::Constant *(*)(llvm::Type *, double)) &llvm::ConstantFP::get,"");
}
NOINLINE void expose_function_155_helper() {
  expose_function(core::magic_name("LLVM-SYS:constantFpGet"),true,(llvm::ConstantFP *(*)(llvm::LLVMContext &, const llvm::APFloat &)) &llvm::ConstantFP::get,"");
}
NOINLINE void expose_function_156_helper() {
  expose_function(core::magic_name("LLVM-SYS:make-global-variable"),true,&llvmo::GlobalVariable_O::make,"(module type is-constant linkage initializer name &optional (insert-before nil) (thread-local-mode 'llvm-sys:not-thread-local))");
}
NOINLINE void expose_function_157_helper() {
  expose_function(core::magic_name("LLVM-SYS:constant-expr-get-in-bounds-get-element-ptr"),true,&llvmo::ConstantExpr_O::getInBoundsGetElementPtr,"");
}
NOINLINE void expose_function_158_helper() {
  expose_function(core::magic_name("LLVM-SYS:block-address-get"),true,&llvmo::BlockAddress_O::get,"(function basic-block)");
}
NOINLINE void expose_function_159_helper() {
  expose_function(core::magic_name("LLVM-SYS:constant-array-get"),true,&llvmo::ConstantArray_O::get,"(type values)");
}
NOINLINE void expose_function_160_helper() {
  expose_function(core::magic_name("LLVM-SYS:constant-data-array-get-uint32"),true,&llvmo::ConstantDataArray_O::getUInt32,"(type values)");
}
NOINLINE void expose_function_161_helper() {
  expose_function(core::magic_name("LLVM-SYS:PassManagerBuilderSetfSizeLevel"),true,&llvmo::PassManagerBuilderSetfSizeLevel,"");
}
NOINLINE void expose_function_162_helper() {
  expose_function(core::magic_name("LLVM-SYS:PassManagerBuilderSetfOptLevel"),true,&llvmo::PassManagerBuilderSetfOptLevel,"");
}
NOINLINE void expose_function_163_helper() {
  expose_function(core::magic_name("LLVM-SYS:PassManagerBuilderSetfInliner"),true,&llvmo::PassManagerBuilderSetfInliner,"");
}
NOINLINE void expose_function_164_helper() {
  expose_function(core::magic_name("make-PassManagerBuilder","LLVM-SYS"),true,&llvmo::PassManagerBuilder_O::make,"");
}
NOINLINE void expose_function_165_helper() {
  expose_function(core::magic_name("make-EngineBuilder","LLVM-SYS"),true,&llvmo::EngineBuilder_O::make,"(module)");
}
NOINLINE void expose_function_166_helper() {
  expose_function(core::magic_name("makePassManager","LLVM-SYS"),true,&llvmo::PassManager_O::make,"");
}
NOINLINE void expose_function_167_helper() {
  expose_function(core::magic_name("makeFunctionPassManager","LLVM-SYS"),true,&llvmo::FunctionPassManager_O::make,"(module)");
}
NOINLINE void expose_function_168_helper() {
  expose_function(core::magic_name("makeTargetLibraryInfo","LLVM-SYS"),true,&llvmo::TargetLibraryInfo_O::make,"(triple)");
}
NOINLINE void expose_function_169_helper() {
  expose_function(core::magic_name("makeDataLayoutPass","LLVM-SYS"),true,&llvmo::DataLayoutPass_O::make,"(module)");
}
NOINLINE void expose_function_170_helper() {
  expose_function(core::magic_name("llvm_sys__module_get_function_list"),true,&llvmo::llvm_sys__module_get_function_list,"");
}
NOINLINE void expose_function_171_helper() {
  expose_function(core::magic_name("make-Module","LLVM-SYS"),true,&llvmo::Module_O::make,"(module-name context)");
}
NOINLINE void expose_function_172_helper() {
  expose_function(core::magic_name("llvm_sys__valid"),true,&llvmo::llvm_sys__valid,"");
}
NOINLINE void expose_function_173_helper() {
  expose_function(core::magic_name("llvm_sys__makeStringGlobal"),true,&llvmo::llvm_sys__makeStringGlobal,"");
}
NOINLINE void expose_function_174_helper() {
  expose_function(core::magic_name("attributeSetGet","LLVM-SYS"),true,(llvm::AttributeSet (*)(llvm::LLVMContext &, unsigned, llvm::ArrayRef<llvm::Attribute::AttrKind>)) & llvm::AttributeSet::get,"");
}
NOINLINE void expose_function_175_helper() {
  expose_function(core::magic_name("llvm_sys__valuep"),true,&llvmo::llvm_sys__valuep,"");
}
NOINLINE void expose_function_176_helper() {
  expose_function(core::magic_name("llvm_sys__parseBitcodeFile"),true,&llvmo::llvm_sys__parseBitcodeFile,"");
}
NOINLINE void expose_function_177_helper() {
  expose_function(core::magic_name("llvm_sys__writeBitcodeToFile"),true,&llvmo::llvm_sys__writeBitcodeToFile,"");
}
NOINLINE void expose_function_178_helper() {
  expose_function(core::magic_name("llvm_sys__verifyFunction"),true,&llvmo::llvm_sys__verifyFunction,"");
}
NOINLINE void expose_function_179_helper() {
  expose_function(core::magic_name("llvm_sys__verifyModule"),true,&llvmo::llvm_sys__verifyModule,"");
}
NOINLINE void expose_function_180_helper() {
  expose_function(core::magic_name("llvm_sys__writeIrToFile"),true,&llvmo::llvm_sys__writeIrToFile,"");
}
NOINLINE void expose_function_181_helper() {
  expose_function(core::magic_name("LLVM-SYS:make-target-options"),true,&llvmo::TargetOptions_O::make,"");
}
NOINLINE void expose_function_182_helper() {
  expose_function(core::magic_name("triple-normalize","LLVM-SYS"),true,&llvm::Triple::normalize,"");
}
NOINLINE void expose_function_183_helper() {
  expose_function(core::magic_name("make-triple","LLVM-SYS"),true,&llvmo::Triple_O::make,"(triple-str)");
}
NOINLINE void expose_function_184_helper() {
  expose_function(core::magic_name("llvm_sys__link_in_module"),true,&llvmo::llvm_sys__link_in_module,"");
}
NOINLINE void expose_function_185_helper() {
  expose_function(core::magic_name("LLVM-SYS:make-linker"),true,&llvmo::Linker_O::make,"");
}
NOINLINE void expose_function_186_helper() {
  expose_function(core::magic_name("LLVM-SYS:get_global_context"),true,&llvmo::LLVMContext_O::get_global_context,"");
}
NOINLINE void expose_function_187_helper() {
  expose_function(core::magic_name("llvm_sys__llvm_value_p"),true,&llvmo::llvm_sys__llvm_value_p,"");
}
NOINLINE void expose_function_188_helper() {
  expose_function(core::magic_name("compiler__setAssociatedFuncs"),true,&llvmo::compiler__setAssociatedFuncs,"");
}
NOINLINE void expose_function_189_helper() {
  expose_function("LLVM-SYS:dw-tag",true,&llvmo::llvm_sys__dwTag,"(tagsym &optional (debug-version llvm-sys:llvmdebug-version8))");
}
NOINLINE void expose_function_190_helper() {
  expose_function("LLVM-SYS:debug-loc-get",true,&llvmo::DebugLoc_O::get,"");
}
NOINLINE void expose_function_191_helper() {
  expose_function(core::magic_name("LLVM-SYS:make-dibuilder"),true,&llvmo::DIBuilder_O::make,"");
}
NOINLINE void expose_function_192_helper() {
  expose_function(core::magic_name("llvm_sys__addGlobalBootFunctionsSizePass"),true,&llvmo::llvm_sys__addGlobalBootFunctionsSizePass,"");
}
NOINLINE void expose_function_193_helper() {
  expose_function(core::magic_name("null-pointer","CFFI-SYS"),true,&cffi::Pointer_O::null_pointer,"");
}
NOINLINE void expose_function_194_helper() {
  expose_function(core::magic_name("make-pointer","CFFI-SYS"),true,&cffi::Pointer_O::make,"");
}
NOINLINE void expose_function_195_helper() {
  expose_function(core::magic_name("cffi_sys__foreign_alloc"),true,&cffi::cffi_sys__foreign_alloc,"");
}
NOINLINE void expose_function_196_helper() {
  expose_function(core::magic_name("cffi_sys__PERCENTforeign_type_size"),true,&cffi::cffi_sys__PERCENTforeign_type_size,"");
}
NOINLINE void expose_function_197_helper() {
  expose_function(core::magic_name("cffi_sys__PERCENTforeign_type_alignment"),true,&cffi::cffi_sys__PERCENTforeign_type_alignment,"");
}
NOINLINE void expose_function_198_helper() {
  expose_function(core::magic_name("cffi_sys__foreign_symbol_pointer"),true,&cffi::cffi_sys__foreign_symbol_pointer,"");
}
NOINLINE void expose_function_199_helper() {
  expose_function(core::magic_name("cffi_sys__PERCENTload_foreign_library"),true,&cffi::cffi_sys__PERCENTload_foreign_library,"");
}
NOINLINE void expose_function_200_helper() {
  expose_function(core::magic_name("cl__get_internal_run_time"),true,&core::cl__get_internal_run_time,"");
}
NOINLINE void expose_function_201_helper() {
  expose_function(core::magic_name("cl__get_internal_real_time"),true,&core::cl__get_internal_real_time,"");
}
NOINLINE void expose_function_202_helper() {
  expose_function(core::magic_name("make-bit-vector-with-fill-ptr","CORE"),true,&core::BitVectorWithFillPtr_O::make,"");
}
NOINLINE void expose_function_203_helper() {
  expose_function(core::magic_name("make-simple-bit-vector","CORE"),true,&core::SimpleBitVector_O::make,"");
}
NOINLINE void expose_function_204_helper() {
  expose_function(core::magic_name("cl__ensure_directories_exist"),true,&core::cl__ensure_directories_exist,"(pathspec)");
}
NOINLINE void expose_function_205_helper() {
  expose_function(core::magic_name("CORE:make-array-objects"),true,&core::ArrayObjects_O::make,"");
}
NOINLINE void expose_function_206_helper() {
  expose_function(core::magic_name("CORE:make-load-time-values"),true,&core::LoadTimeValues_O::make,"");
}
NOINLINE void expose_function_207_helper() {
  expose_function(core::magic_name("core__load_time_values_dump_symbols"),true,&core::core__load_time_values_dump_symbols,"(name-or-ltv &optional indices)");
}
NOINLINE void expose_function_208_helper() {
  expose_function(core::magic_name("core__load_time_values_dump_values"),true,&core::core__load_time_values_dump_values,"(name-or-ltv &optional indices)");
}
NOINLINE void expose_function_209_helper() {
  expose_function(core::magic_name("core__load_time_values_ids"),true,&core::core__load_time_values_ids,"");
}
NOINLINE void expose_function_210_helper() {
  expose_function(core::magic_name("core__lookup_load_time_symbol"),true,&core::core__lookup_load_time_symbol,"(name idx)");
}
NOINLINE void expose_function_211_helper() {
  expose_function(core::magic_name("core__lookup_load_time_value"),true,&core::core__lookup_load_time_value,"(name idx)");
}
NOINLINE void expose_function_212_helper() {
  expose_function(core::magic_name("core__load_time_value_array"),true,&core::core__load_time_value_array,"(name &optional (data-size 0) (symbol-size 0))");
}
NOINLINE void expose_function_213_helper() {
  expose_function(core::magic_name("core__set_run_time_values_vector"),true,&core::core__set_run_time_values_vector,"(name)");
}
NOINLINE void expose_function_214_helper() {
  expose_function(core::magic_name("cl__vector_push_extend"),true,&core::cl__vector_push_extend,"(newElement vector &optional (exension 16))");
}
NOINLINE void expose_function_215_helper() {
  expose_function(core::magic_name("cl__vector_push"),true,&core::cl__vector_push,"(newElement vector)");
}
NOINLINE void expose_function_216_helper() {
  expose_function(core::magic_name("core__adjust_vector"),true,&core::core__adjust_vector,"(array dimensions initial-element initial-contents)");
}
NOINLINE void expose_function_217_helper() {
  expose_function(core::magic_name("core__make_vector"),true,&core::core__make_vector,"(element-type dimension &optional adjustable (fill-pointer t) displaced-to displaced-index-offset initial-element initial-contents)");
}
NOINLINE void expose_function_218_helper() {
  expose_function(core::magic_name("cl__vector"),true,&core::cl__vector,"(&rest args)");
}
NOINLINE void expose_function_219_helper() {
  expose_function(core::magic_name("core__setf_subseq"),true,&core::core__setf_subseq,"(sequence start end subseq)");
}
NOINLINE void expose_function_220_helper() {
  expose_function(core::magic_name("core__elt_set"),true,&core::core__elt_set,"(seq index val)");
}
NOINLINE void expose_function_221_helper() {
  expose_function(core::magic_name("core__setf_elt"),true,&core::core__setf_elt,"(sequence index value)");
}
NOINLINE void expose_function_222_helper() {
  expose_function(core::magic_name("cl__copy_seq"),true,&core::cl__copy_seq,"(seq)");
}
NOINLINE void expose_function_223_helper() {
  expose_function(core::magic_name("cl__subseq"),true,&core::cl__subseq,"(sequence start &optional end)");
}
NOINLINE void expose_function_224_helper() {
  expose_function(core::magic_name("cl__nreverse"),true,&core::cl__nreverse,"(seq)");
}
NOINLINE void expose_function_225_helper() {
  expose_function(core::magic_name("cl__reverse"),true,&core::cl__reverse,"(seq)");
}
NOINLINE void expose_function_226_helper() {
  expose_function(core::magic_name("cl__make_sequence"),true,&core::cl__make_sequence,"(type size &key (initial-element nil iesp))");
}
NOINLINE void expose_function_227_helper() {
  expose_function(core::magic_name("cl__elt"),true,&core::cl__elt,"(sequence index)");
}
NOINLINE void expose_function_228_helper() {
  expose_function(core::magic_name("cl__length"),true,&core::cl__length,"(arg)");
}
NOINLINE void expose_function_229_helper() {
  expose_function(core::magic_name("coerce-to-package","CORE"),true,&core::coerce_to_package,"");
}
NOINLINE void expose_function_230_helper() {
  expose_function(core::magic_name("path-designator","CORE"),true,&core::core__path_designator,"");
}
NOINLINE void expose_function_231_helper() {
  expose_function(core::magic_name("cl__parse_integer"),true,&core::cl__parse_integer,"(string &key (start 0) end (radix 10) junk-allowed)");
}
NOINLINE void expose_function_232_helper() {
  expose_function(core::magic_name("core__schar_set"),true,&core::core__schar_set,"(str index c)");
}
NOINLINE void expose_function_233_helper() {
  expose_function(core::magic_name("core__char_set"),true,&core::core__char_set,"(str index c)");
}
NOINLINE void expose_function_234_helper() {
  expose_function(core::magic_name("cl__schar"),true,&core::cl__schar,"(str index)");
}
NOINLINE void expose_function_235_helper() {
  expose_function(core::magic_name("cl__make_string"),true,&core::cl__make_string,"(size &key initial-element (element-type 'character))");
}
NOINLINE void expose_function_236_helper() {
  expose_function(core::magic_name("cl__string_not_lessp"),true,&core::cl__string_not_lessp,"(strdes1 strdes2 &key (start1 0) end1 (start2 0) end2)");
}
NOINLINE void expose_function_237_helper() {
  expose_function(core::magic_name("cl__string_not_greaterp"),true,&core::cl__string_not_greaterp,"(strdes1 strdes2 &key (start1 0) end1 (start2 0) end2)");
}
NOINLINE void expose_function_238_helper() {
  expose_function(core::magic_name("cl__string_greaterp"),true,&core::cl__string_greaterp,"(strdes1 strdes2 &key (start1 0) end1 (start2 0) end2)");
}
NOINLINE void expose_function_239_helper() {
  expose_function(core::magic_name("cl__string_lessp"),true,&core::cl__string_lessp,"(strdes1 strdes2 &key (start1 0) end1 (start2 0) end2)");
}
NOINLINE void expose_function_240_helper() {
  expose_function(core::magic_name("cl__string_not_equal"),true,&core::cl__string_not_equal,"(strdes1 strdes2 &key (start1 0) end1 (start2 0) end2)");
}
NOINLINE void expose_function_241_helper() {
  expose_function(core::magic_name("cl__string_equal"),true,&core::cl__string_equal,"(strdes1 strdes2 &key (start1 0) end1 (start2 0) end2)");
}
NOINLINE void expose_function_242_helper() {
  expose_function(core::magic_name("cl__string_GE_"),true,&core::cl__string_GE_,"(strdes1 strdes2 &key (start1 0) end1 (start2 0) end2)");
}
NOINLINE void expose_function_243_helper() {
  expose_function(core::magic_name("cl__string_LE_"),true,&core::cl__string_LE_,"(strdes1 strdes2 &key (start1 0) end1 (start2 0) end2)");
}
NOINLINE void expose_function_244_helper() {
  expose_function(core::magic_name("cl__string_GT_"),true,&core::cl__string_GT_,"(strdes1 strdes2 &key (start1 0) end1 (start2 0) end2)");
}
NOINLINE void expose_function_245_helper() {
  expose_function(core::magic_name("cl__string_LT_"),true,&core::cl__string_LT_,"(strdes1 strdes2 &key (start1 0) end1 (start2 0) end2)");
}
NOINLINE void expose_function_246_helper() {
  expose_function(core::magic_name("cl__string_NE_"),true,&core::cl__string_NE_,"(strdes1 strdes2 &key (start1 0) end1 (start2 0) end2)");
}
NOINLINE void expose_function_247_helper() {
  expose_function(core::magic_name("cl__string_EQ_"),true,&core::cl__string_EQ_,"(strdes1 strdes2 &key (start1 0) end1 (start2 0) end2)");
}
NOINLINE void expose_function_248_helper() {
  expose_function(core::magic_name("cl__string"),true,&core::cl__string,"(arg)");
}
NOINLINE void expose_function_249_helper() {
  expose_function(core::magic_name("cl__nstring_downcase"),true,&core::cl__nstring_downcase,"(arg)");
}
NOINLINE void expose_function_250_helper() {
  expose_function(core::magic_name("cl__nstring_upcase"),true,&core::cl__nstring_upcase,"(arg)");
}
NOINLINE void expose_function_251_helper() {
  expose_function(core::magic_name("cl__string_downcase"),true,&core::cl__string_downcase,"(arg)");
}
NOINLINE void expose_function_252_helper() {
  expose_function(core::magic_name("cl__char"),true,&core::cl__char,"(str idx)");
}
NOINLINE void expose_function_253_helper() {
  expose_function(core::magic_name("cl__string_upcase"),true,&core::cl__string_upcase,"(arg)");
}
NOINLINE void expose_function_254_helper() {
  expose_function(core::magic_name("cl__string_right_trim"),true,&core::cl__string_right_trim,"(charbag str)");
}
NOINLINE void expose_function_255_helper() {
  expose_function(core::magic_name("cl__string_left_trim"),true,&core::cl__string_left_trim,"(charbag str)");
}
NOINLINE void expose_function_256_helper() {
  expose_function(core::magic_name("cl__string_trim"),true,&core::cl__string_trim,"(charbag str)");
}
NOINLINE void expose_function_257_helper() {
  expose_function(core::magic_name("cl__code_char"),true,&core::cl__code_char,"(och)");
}
NOINLINE void expose_function_258_helper() {
  expose_function(core::magic_name("cl__char_int"),true,&core::cl__char_int,"(och)");
}
NOINLINE void expose_function_259_helper() {
  expose_function(core::magic_name("cl__char_code"),true,&core::cl__char_code,"(och)");
}
NOINLINE void expose_function_260_helper() {
  expose_function(core::magic_name("cl__char_name"),true,&core::cl__char_name,"(och)");
}
NOINLINE void expose_function_261_helper() {
  expose_function(core::magic_name("cl__name_char"),true,&core::cl__name_char,"(sname)");
}
NOINLINE void expose_function_262_helper() {
  expose_function(core::magic_name("cl__digit_char_p"),true,&core::cl__digit_char_p,"(c &optional (radix 10))");
}
NOINLINE void expose_function_263_helper() {
  expose_function(core::magic_name("cl__alpha_char_p"),true,&core::cl__alpha_char_p,"(ch)");
}
NOINLINE void expose_function_264_helper() {
  expose_function(core::magic_name("cl__standard_char_p"),true,&core::cl__standard_char_p,"(ch)");
}
NOINLINE void expose_function_265_helper() {
  expose_function(core::magic_name("cl__char_equal"),true,&core::cl__char_equal,"(&rest args)");
}
NOINLINE void expose_function_266_helper() {
  expose_function(core::magic_name("cl__char_not_equal"),true,&core::cl__char_not_equal,"(&rest args)");
}
NOINLINE void expose_function_267_helper() {
  expose_function(core::magic_name("cl__char_EQ_"),true,&core::cl__char_EQ_,"(&rest args)");
}
NOINLINE void expose_function_268_helper() {
  expose_function(core::magic_name("cl__char_NE_"),true,&core::cl__char_NE_,"(&rest args)");
}
NOINLINE void expose_function_269_helper() {
  expose_function(core::magic_name("cl__char_not_lessp"),true,&core::cl__char_not_lessp,"(&rest args)");
}
NOINLINE void expose_function_270_helper() {
  expose_function(core::magic_name("cl__char_not_greaterp"),true,&core::cl__char_not_greaterp,"(&rest args)");
}
NOINLINE void expose_function_271_helper() {
  expose_function(core::magic_name("cl__char_greaterp"),true,&core::cl__char_greaterp,"(&rest args)");
}
NOINLINE void expose_function_272_helper() {
  expose_function(core::magic_name("cl__char_lessp"),true,&core::cl__char_lessp,"(&rest args)");
}
NOINLINE void expose_function_273_helper() {
  expose_function(core::magic_name("cl__char_GE_"),true,&core::cl__char_GE_,"(&rest args)");
}
NOINLINE void expose_function_274_helper() {
  expose_function(core::magic_name("cl__char_LE_"),true,&core::cl__char_LE_,"(&rest args)");
}
NOINLINE void expose_function_275_helper() {
  expose_function(core::magic_name("cl__char_GT_"),true,&core::cl__char_GT_,"(&rest args)");
}
NOINLINE void expose_function_276_helper() {
  expose_function(core::magic_name("cl__char_LT_"),true,&core::cl__char_LT_,"(&rest args)");
}
NOINLINE void expose_function_277_helper() {
  expose_function(core::magic_name("cl__char_downcase"),true,&core::cl__char_downcase,"(char)");
}
NOINLINE void expose_function_278_helper() {
  expose_function(core::magic_name("cl__char_upcase"),true,&core::cl__char_upcase,"(char)");
}
NOINLINE void expose_function_279_helper() {
  expose_function(core::magic_name("cl__alphanumericp"),true,&core::cl__alphanumericp,"(char)");
}
NOINLINE void expose_function_280_helper() {
  expose_function(core::magic_name("cl__both_case_p"),true,&core::cl__both_case_p,"(arg)");
}
NOINLINE void expose_function_281_helper() {
  expose_function(core::magic_name("cl__upper_case_p"),true,&core::cl__upper_case_p,"(arg)");
}
NOINLINE void expose_function_282_helper() {
  expose_function(core::magic_name("cl__lower_case_p"),true,&core::cl__lower_case_p,"(arg)");
}
NOINLINE void expose_function_283_helper() {
  expose_function(core::magic_name("core__structure_subtypep"),true,&core::core__structure_subtypep,"(x y)");
}
NOINLINE void expose_function_284_helper() {
  expose_function(core::magic_name("core__structurep"),true,&core::core__structurep,"(arg)");
}
NOINLINE void expose_function_285_helper() {
  expose_function(core::magic_name("core__structure_set"),true,&core::core__structure_set,"(struct type idx val)");
}
NOINLINE void expose_function_286_helper() {
  expose_function(core::magic_name("core__structure_ref"),true,&core::core__structure_ref,"(obj name idx)");
}
NOINLINE void expose_function_287_helper() {
  expose_function(core::magic_name("cl__copy_structure"),true,&core::cl__copy_structure,"(arg)");
}
NOINLINE void expose_function_288_helper() {
  expose_function(core::magic_name("core__make_structure"),true,&core::core__make_structure,"(type &rest slot_values)");
}
NOINLINE void expose_function_289_helper() {
  expose_function(core::magic_name("core__ensure_single_dispatch_method"),true,&core::core__ensure_single_dispatch_method,"(gfname receiver-class &key lambda-list-handler declares (docstring \"\") body )");
}
NOINLINE void expose_function_290_helper() {
  expose_function(core::magic_name("core__ensure_single_dispatch_generic_function"),true,&core::core__ensure_single_dispatch_generic_function,"(gfname llhandler)");
}
NOINLINE void expose_function_291_helper() {
  expose_function(core::magic_name("core__nread"),true,&core::core__nread,"(sin &optional (eof-error-p t) eof-value)");
}
NOINLINE void expose_function_292_helper() {
  expose_function(core::magic_name("core__ensure_documentation"),true,&core::core__ensure_documentation,"(sub-key symbol value)");
}
NOINLINE void expose_function_293_helper() {
  expose_function(core::magic_name("ext__annotate"),true,&core::ext__annotate,"(object key sub-key value)");
}
NOINLINE void expose_function_294_helper() {
  expose_function(core::magic_name("core__rem_record_field"),true,&core::core__rem_record_field,"(record key sub-key)");
}
NOINLINE void expose_function_295_helper() {
  expose_function(core::magic_name("core__set_record_field"),true,&core::core__set_record_field,"(record key sub-key value)");
}
NOINLINE void expose_function_296_helper() {
  expose_function(core::magic_name("core__record_field"),true,&core::core__record_field,"(record key sub-key)");
}
NOINLINE void expose_function_297_helper() {
  expose_function(core::magic_name("core__record_cons"),true,&core::core__record_cons,"(record key sub-key)");
}
NOINLINE void expose_function_298_helper() {
  expose_function(core::magic_name("core__backquote_remove_tokens"),true,&core::core__backquote_remove_tokens,"(x)");
}
NOINLINE void expose_function_299_helper() {
  expose_function(core::magic_name("core__backquote_attach_conses"),true,&core::core__backquote_attach_conses,"(items result)");
}
NOINLINE void expose_function_300_helper() {
  expose_function(core::magic_name("core__backquote_attach_append"),true,&core::core__backquote_attach_append,"(op item result)");
}
NOINLINE void expose_function_301_helper() {
  expose_function(core::magic_name("core__backquote_null_or_quoted"),true,&core::core__backquote_null_or_quoted,"(x)");
}
NOINLINE void expose_function_302_helper() {
  expose_function(core::magic_name("core__backquote_simplify_args"),true,&core::core__backquote_simplify_args,"(x)");
}
NOINLINE void expose_function_303_helper() {
  expose_function(core::magic_name("core__backquote_simplify"),true,&core::core__backquote_simplify,"(x)");
}
NOINLINE void expose_function_304_helper() {
  expose_function(core::magic_name("core__backquote_maptree"),true,&core::core__backquote_maptree,"(op x)");
}
NOINLINE void expose_function_305_helper() {
  expose_function(core::magic_name("core__backquote_frob"),true,&core::core__backquote_frob,"(x)");
}
NOINLINE void expose_function_306_helper() {
  expose_function(core::magic_name("core__backquote_splicing_frob"),true,&core::core__backquote_splicing_frob,"(x)");
}
NOINLINE void expose_function_307_helper() {
  expose_function(core::magic_name("core__backquote_bracket"),true,&core::core__backquote_bracket,"(x)");
}
NOINLINE void expose_function_308_helper() {
  expose_function(core::magic_name("core__backquote_process"),true,&core::core__backquote_process,"(ox)");
}
NOINLINE void expose_function_309_helper() {
  expose_function(core::magic_name("core__backquote_completely_process"),true,&core::core__backquote_completely_process,"(x)");
}
NOINLINE void expose_function_310_helper() {
  expose_function(core::magic_name("core__backquote_append"),true,&core::core__backquote_append,"(&rest lists)");
}
NOINLINE void expose_function_311_helper() {
  expose_function(core::magic_name("cl__format"),true,&core::cl__format,"(destination control &rest args)");
}
NOINLINE void expose_function_312_helper() {
  expose_function(core::magic_name("core__bformat"),true,&core::core__bformat,"(destination control &rest args)");
}
NOINLINE void expose_function_313_helper() {
  expose_function(core::magic_name("core__rem_sysprop"),true,&core::core__rem_sysprop,"(key area)");
}
NOINLINE void expose_function_314_helper() {
  expose_function(core::magic_name("core__get_sysprop"),true,&core::core__get_sysprop,"(key area)");
}
NOINLINE void expose_function_315_helper() {
  expose_function(core::magic_name("core__put_sysprop"),true,&core::core__put_sysprop,"(key area value)");
}
NOINLINE void expose_function_316_helper() {
  expose_function(core::magic_name("core__iterator_step"),true,&core::core__iterator_step,"(it)");
}
NOINLINE void expose_function_317_helper() {
  expose_function(core::magic_name("core__iterator_unsafe_element"),true,&core::core__iterator_unsafe_element,"(it)");
}
NOINLINE void expose_function_318_helper() {
  expose_function(core::magic_name("CORE:make-array-displaced"),true,&core::ArrayDisplaced_O::make,"");
}
NOINLINE void expose_function_319_helper() {
  expose_function(core::magic_name("core__condition_message"),true,&core::core__condition_message,"(c)");
}
NOINLINE void expose_function_320_helper() {
  expose_function(core::magic_name("cl__make_condition"),true,&core::cl__make_condition,"(type &rest slot-initializations)");
}
NOINLINE void expose_function_321_helper() {
  expose_function(core::magic_name("cl__make_symbol"),true,&core::cl__make_symbol,"(name)");
}
NOINLINE void expose_function_322_helper() {
  expose_function(core::magic_name("core__symbol_value_address"),true,&core::core__symbol_value_address,"(arg)");
}
NOINLINE void expose_function_323_helper() {
  expose_function(core::magic_name("cl__symbol_value"),true,&core::cl__symbol_value,"(arg)");
}
NOINLINE void expose_function_324_helper() {
  expose_function(core::magic_name("cl__symbol_name"),true,&core::cl__symbol_name,"(arg)");
}
NOINLINE void expose_function_325_helper() {
  expose_function(core::magic_name("cl__symbol_function"),true,&core::cl__symbol_function,"(arg)");
}
NOINLINE void expose_function_326_helper() {
  expose_function(core::magic_name("cl__symbol_package"),true,&core::cl__symbol_package,"(arg)");
}
NOINLINE void expose_function_327_helper() {
  expose_function(core::magic_name("cl__boundp"),true,&core::cl__boundp,"(arg)");
}
NOINLINE void expose_function_328_helper() {
  expose_function(core::magic_name("core__putprop"),true,&core::core__putprop,"(sym val indicator)");
}
NOINLINE void expose_function_329_helper() {
  expose_function(core::magic_name("core__setf_symbol_plist"),true,&core::core__setf_symbol_plist,"(sym plist)");
}
NOINLINE void expose_function_330_helper() {
  expose_function(core::magic_name("cl__get"),true,&core::cl__get,"(sym indicator &optional default)");
}
NOINLINE void expose_function_331_helper() {
  expose_function(core::magic_name("cl__symbol_plist"),true,&core::cl__symbol_plist,"(sym)");
}
NOINLINE void expose_function_332_helper() {
  expose_function(core::magic_name("core__alist_asString"),true,&core::core__alist_asString,"");
}
NOINLINE void expose_function_333_helper() {
  expose_function(core::magic_name("core__alist_get"),true,&core::core__alist_get,"");
}
NOINLINE void expose_function_334_helper() {
  expose_function(core::magic_name("core__alist_push"),true,&core::core__alist_push,"");
}
NOINLINE void expose_function_335_helper() {
  expose_function(core::magic_name("core__alist_erase"),true,&core::core__alist_erase,"");
}
NOINLINE void expose_function_336_helper() {
  expose_function(core::magic_name("core__append2"),true,&core::core__append2,"(l1 l2)");
}
NOINLINE void expose_function_337_helper() {
  expose_function(core::magic_name("cl__make_list"),true,&core::cl__make_list,"(osize &key initial_element)");
}
NOINLINE void expose_function_338_helper() {
  expose_function(core::magic_name("cl__rplacd"),true,&core::cl__rplacd,"(c o)");
}
NOINLINE void expose_function_339_helper() {
  expose_function(core::magic_name("cl__rplaca"),true,&core::cl__rplaca,"(c o)");
}
NOINLINE void expose_function_340_helper() {
  expose_function(core::magic_name("cl__cons"),true,&core::cl__cons,"(object1 object2)");
}
NOINLINE void expose_function_341_helper() {
  expose_function(core::magic_name("core__rem_f"),true,&core::core__rem_f,"(plist indicator)");
}
NOINLINE void expose_function_342_helper() {
  expose_function(core::magic_name("cl__getf"),true,&core::cl__getf,"(plist indicator &optional default-value)");
}
NOINLINE void expose_function_343_helper() {
  expose_function(core::magic_name("core__put_f"),true,&core::core__put_f,"(plist value indicator)");
}
NOINLINE void expose_function_344_helper() {
  expose_function(core::magic_name("core__aset"),true,&core::core__aset,"(array &rest indices-value)");
}
NOINLINE void expose_function_345_helper() {
  expose_function(core::magic_name("core__copy_subarray"),true,&core::core__copy_subarray,"(out outStart in inStart len)");
}
NOINLINE void expose_function_346_helper() {
  expose_function(core::magic_name("cl__upgraded_array_element_type"),true,&core::cl__upgraded_array_element_type,"(core::type &optional core::env)");
}
NOINLINE void expose_function_347_helper() {
  expose_function(core::magic_name("cl__array_displacement"),true,&core::cl__array_displacement,"(core::array)");
}
NOINLINE void expose_function_348_helper() {
  expose_function(core::magic_name("cl__fill_pointer"),true,&core::cl__fill_pointer,"(vec)");
}
NOINLINE void expose_function_349_helper() {
  expose_function(core::magic_name("core__calculate_missing_common_lisp_symbols"),true,&core::core__calculate_missing_common_lisp_symbols,"");
}
NOINLINE void expose_function_350_helper() {
  expose_function(core::magic_name("cl__package_name"),true,&core::cl__package_name,"(pkg)");
}
NOINLINE void expose_function_351_helper() {
  expose_function(core::magic_name("cl__package_used_by_list"),true,&core::cl__package_used_by_list,"(pkg)");
}
NOINLINE void expose_function_352_helper() {
  expose_function(core::magic_name("cl__package_use_list"),true,&core::cl__package_use_list,"(package-designator)");
}
NOINLINE void expose_function_353_helper() {
  expose_function(core::magic_name("cl__gentemp"),true,&core::cl__gentemp,"(&optional prefix (package *package*))");
}
NOINLINE void expose_function_354_helper() {
  expose_function(core::magic_name("cl__shadowing_import"),true,&core::cl__shadowing_import,"(symbol-names-desig &optional (package-desig *package*))");
}
NOINLINE void expose_function_355_helper() {
  expose_function(core::magic_name("cl__shadow"),true,&core::cl__shadow,"(symbol-names-desig &optional (package_desig *package*))");
}
NOINLINE void expose_function_356_helper() {
  expose_function(core::magic_name("cl__import"),true,&core::cl__import,"(symbols-desig &optional (package-desig *package*))");
}
NOINLINE void expose_function_357_helper() {
  expose_function(core::magic_name("cl__package_shadowing_symbols"),true,&core::cl__package_shadowing_symbols,"(package_desig)");
}
NOINLINE void expose_function_358_helper() {
  expose_function(core::magic_name("cl__delete_package"),true,&core::cl__delete_package,"(packages-to-unuse-desig &optional (package-desig *package*))");
}
NOINLINE void expose_function_359_helper() {
  expose_function(core::magic_name("cl__unuse_package"),true,&core::cl__unuse_package,"(packages-to-unuse-desig &optional (package-desig *package*))");
}
NOINLINE void expose_function_360_helper() {
  expose_function(core::magic_name("cl__use_package"),true,&core::cl__use_package,"(packages-to-use-desig &optional (package-desig *package*))");
}
NOINLINE void expose_function_361_helper() {
  expose_function(core::magic_name("cl__list_all_packages"),true,&core::cl__list_all_packages,"");
}
NOINLINE void expose_function_362_helper() {
  expose_function(core::magic_name("cl__make_package"),true,&core::cl__make_package,"(package-name &key nicknames (use (list \"CL\")))");
}
NOINLINE void expose_function_363_helper() {
  expose_function(core::magic_name("cl__find_symbol"),true,&core::cl__find_symbol,"(sym &optional (package *package*))");
}
NOINLINE void expose_function_364_helper() {
  expose_function(core::magic_name("cl__unintern"),true,&core::cl__unintern,"(symbol &optional (package *package*))");
}
NOINLINE void expose_function_365_helper() {
  expose_function(core::magic_name("cl__package_nicknames"),true,&core::cl__package_nicknames,"(pkg)");
}
NOINLINE void expose_function_366_helper() {
  expose_function(core::magic_name("cl__rename_package"),true,&core::cl__rename_package,"(package new-name &optional nick-names)");
}
NOINLINE void expose_function_367_helper() {
  expose_function(core::magic_name("CORE:base_string_concatenate"),true,&core::core__base_string_concatenate,"(&va-rest args)");
}
NOINLINE void expose_function_368_helper() {
  expose_function(core::magic_name("core__search_string"),true,&core::core__search_string,"(str1 start1 end1 str2 start2 end2)");
}
NOINLINE void expose_function_369_helper() {
  expose_function(core::magic_name("core__logical_pathname_p"),true,&core::core__logical_pathname_p,"(arg)");
}
NOINLINE void expose_function_370_helper() {
  expose_function(core::magic_name("cl__pathnamep"),true,&core::cl__pathnamep,"(arg)");
}
NOINLINE void expose_function_371_helper() {
  expose_function(core::magic_name("core__proper_list_p"),true,&core::core__proper_list_p,"(arg)");
}
NOINLINE void expose_function_372_helper() {
  expose_function(core::magic_name("core__generic_function_p"),true,&core::core__generic_function_p,"(arg)");
}
NOINLINE void expose_function_373_helper() {
  expose_function(core::magic_name("cl__compiled_function_p"),true,&core::cl__compiled_function_p,"(arg)");
}
NOINLINE void expose_function_374_helper() {
  expose_function(core::magic_name("cl__simple_vector_p"),true,&core::cl__simple_vector_p,"(arg)");
}
NOINLINE void expose_function_375_helper() {
  expose_function(core::magic_name("cl__simple_bit_vector_p"),true,&core::cl__simple_bit_vector_p,"(arg)");
}
NOINLINE void expose_function_376_helper() {
  expose_function(core::magic_name("core__external_object_p"),true,&core::core__external_object_p,"(arg)");
}
NOINLINE void expose_function_377_helper() {
  expose_function(core::magic_name("core__activation_frame_p"),true,&core::core__activation_frame_p,"(arg)");
}
NOINLINE void expose_function_378_helper() {
  expose_function(core::magic_name("core__single_dispatch_generic_function_p"),true,&core::core__single_dispatch_generic_function_p,"(arg)");
}
NOINLINE void expose_function_379_helper() {
  expose_function(core::magic_name("core__arrayp"),true,&core::core__arrayp,"");
}
NOINLINE void expose_function_380_helper() {
  expose_function(core::magic_name("core__structure_object_p"),true,&core::core__structure_object_p,"(arg)");
}
NOINLINE void expose_function_381_helper() {
  expose_function(core::magic_name("cl__readtablep"),true,&core::cl__readtablep,"(arg)");
}
NOINLINE void expose_function_382_helper() {
  expose_function(core::magic_name("cl__hash_table_p"),true,&core::cl__hash_table_p,"(arg)");
}
NOINLINE void expose_function_383_helper() {
  expose_function(core::magic_name("cl__bit_vector_p"),true,&core::cl__bit_vector_p,"(arg)");
}
NOINLINE void expose_function_384_helper() {
  expose_function(core::magic_name("core__path_p"),true,&core::core__path_p,"(arg)");
}
NOINLINE void expose_function_385_helper() {
  expose_function(core::magic_name("cl__keywordp"),true,&core::cl__keywordp,"");
}
NOINLINE void expose_function_386_helper() {
  expose_function(core::magic_name("cl__integerp"),true,&core::cl__integerp,"");
}
NOINLINE void expose_function_387_helper() {
  expose_function(core::magic_name("cl__vectorp"),true,&core::cl__vectorp,"");
}
NOINLINE void expose_function_388_helper() {
  expose_function(core::magic_name("cl__floatp"),true,&core::cl__floatp,"");
}
NOINLINE void expose_function_389_helper() {
  expose_function(core::magic_name("cl__realp"),true,&core::cl__realp,"");
}
NOINLINE void expose_function_390_helper() {
  expose_function(core::magic_name("core__single_float_p"),true,&core::core__single_float_p,"(arg)");
}
NOINLINE void expose_function_391_helper() {
  expose_function(core::magic_name("core__short_float_p"),true,&core::core__short_float_p,"(arg)");
}
NOINLINE void expose_function_392_helper() {
  expose_function(core::magic_name("core__long_float_p"),true,&core::core__long_float_p,"(arg)");
}
NOINLINE void expose_function_393_helper() {
  expose_function(core::magic_name("cl__rationalp"),true,&core::cl__rationalp,"");
}
NOINLINE void expose_function_394_helper() {
  expose_function(core::magic_name("cl__random_state_p"),true,&core::cl__random_state_p,"(arg)");
}
NOINLINE void expose_function_395_helper() {
  expose_function(core::magic_name("core__ratio_p"),true,&core::core__ratio_p,"(arg)");
}
NOINLINE void expose_function_396_helper() {
  expose_function(core::magic_name("cl__complexp"),true,&core::cl__complexp,"");
}
NOINLINE void expose_function_397_helper() {
  expose_function(core::magic_name("cl__numberp"),true,&core::cl__numberp,"(arg)");
}
NOINLINE void expose_function_398_helper() {
  expose_function(core::magic_name("core__lambda_list_handler_p"),true,&core::core__lambda_list_handler_p,"(arg)");
}
NOINLINE void expose_function_399_helper() {
  expose_function(core::magic_name("core__cxx_instance_p"),true,&core::core__cxx_instance_p,"(arg)");
}
NOINLINE void expose_function_400_helper() {
  expose_function(core::magic_name("core__built_in_class_p"),true,&core::core__built_in_class_p,"(arg)");
}
NOINLINE void expose_function_401_helper() {
  expose_function(core::magic_name("clos__classp"),true,&core::clos__classp,"");
}
NOINLINE void expose_function_402_helper() {
  expose_function(core::magic_name("cl__packagep"),true,&core::cl__packagep,"(arg)");
}
NOINLINE void expose_function_403_helper() {
  expose_function(core::magic_name("core__interpreted_function_p"),true,&core::core__interpreted_function_p,"(arg)");
}
NOINLINE void expose_function_404_helper() {
  expose_function(core::magic_name("cl__functionp"),true,&core::cl__functionp,"(arg)");
}
NOINLINE void expose_function_405_helper() {
  expose_function(core::magic_name("core__double_float_p"),true,&core::core__double_float_p,"(arg)");
}
NOINLINE void expose_function_406_helper() {
  expose_function(core::magic_name("core__simple_string_p"),true,&core::core__simple_string_p,"(arg)");
}
NOINLINE void expose_function_407_helper() {
  expose_function(core::magic_name("cl__stringp"),true,&core::cl__stringp,"");
}
NOINLINE void expose_function_408_helper() {
  expose_function(core::magic_name("core__bignump"),true,&core::core__bignump,"");
}
NOINLINE void expose_function_409_helper() {
  expose_function(core::magic_name("core__base_char_p"),true,&core::core__base_char_p,"(arg)");
}
NOINLINE void expose_function_410_helper() {
  expose_function(core::magic_name("core__test_val"),true,&core::core__test_val,"(arg)");
}
NOINLINE void expose_function_411_helper() {
  expose_function(core::magic_name("core__is_array"),true,&core::core__is_array,"(arg)");
}
NOINLINE void expose_function_412_helper() {
  expose_function(core::magic_name("core__is_string"),true,&core::core__is_string,"(arg)");
}
NOINLINE void expose_function_413_helper() {
  expose_function(core::magic_name("core__subclassp"),true,&core::core__subclassp,"(low high)");
}
NOINLINE void expose_function_414_helper() {
  expose_function(core::magic_name("core__allocate_raw_class"),true,&core::core__allocate_raw_class,"(original meta-class slots &optional name)");
}
NOINLINE void expose_function_415_helper() {
  expose_function(core::magic_name("core__inherit_default_allocator"),true,&core::core__inherit_default_allocator,"(class directSuperclasses)");
}
NOINLINE void expose_function_416_helper() {
  expose_function(core::magic_name("cl__reader_error"),true,&core::cl__reader_error,"(sourceFileName lineno functionName fmt fmtargs stream)");
}
NOINLINE void expose_function_417_helper() {
  expose_function(core::magic_name("core__wrong_index"),true,&core::core__wrong_index,"(source-file lineno function narg value type)");
}
NOINLINE void expose_function_418_helper() {
  expose_function(core::magic_name("core__wrong_type_nth_arg"),true,&core::core__wrong_type_nth_arg,"(source-file lineno function narg value type)");
}
NOINLINE void expose_function_419_helper() {
  expose_function(core::magic_name("core__wrong_type_argument"),true,&core::core__wrong_type_argument,"(source-file lineno function narg value type)");
}
NOINLINE void expose_function_420_helper() {
  expose_function(core::magic_name("core__signal_simple_error"),true,&core::core__signal_simple_error,"(base-condition continue-message format-control format-args &rest args)");
}
NOINLINE void expose_function_421_helper() {
  expose_function(core::magic_name("cl__translate_logical_pathname"),true,&core::cl__translate_logical_pathname,"(source &key)");
}
NOINLINE void expose_function_422_helper() {
  expose_function(core::magic_name("cl__translate_pathname"),true,&core::cl__translate_pathname,"(source from to &key ((:case scase) :local))");
}
NOINLINE void expose_function_423_helper() {
  expose_function(core::magic_name("core__pathname_translations"),true,&core::core__pathname_translations,"(&optional (host nil hostp) set)");
}
NOINLINE void expose_function_424_helper() {
  expose_function(core::magic_name("cl__pathname_match_p"),true,&core::cl__pathname_match_p,"(tpath tmask)");
}
NOINLINE void expose_function_425_helper() {
  expose_function(core::magic_name("cl__enough_namestring"),true,&core::cl__enough_namestring,"(tpath &optional (defaults (core::safe-default-pathname-defaults)))");
}
NOINLINE void expose_function_426_helper() {
  expose_function(core::magic_name("cl__host_namestring"),true,&core::cl__host_namestring,"(tpname)");
}
NOINLINE void expose_function_427_helper() {
  expose_function(core::magic_name("cl__directory_namestring"),true,&core::cl__directory_namestring,"(tpname)");
}
NOINLINE void expose_function_428_helper() {
  expose_function(core::magic_name("cl__file_namestring"),true,&core::cl__file_namestring,"(tpname)");
}
NOINLINE void expose_function_429_helper() {
  expose_function(core::magic_name("cl__pathname_version"),true,&core::cl__pathname_version,"(pname)");
}
NOINLINE void expose_function_430_helper() {
  expose_function(core::magic_name("cl__pathname_type"),true,&core::cl__pathname_type,"(pname &key ((:case scase) :local))");
}
NOINLINE void expose_function_431_helper() {
  expose_function(core::magic_name("cl__pathname_name"),true,&core::cl__pathname_name,"(pname &key ((:case scase) :local))");
}
NOINLINE void expose_function_432_helper() {
  expose_function(core::magic_name("cl__pathname_directory"),true,&core::cl__pathname_directory,"(pname &key ((:case scase) :local))");
}
NOINLINE void expose_function_433_helper() {
  expose_function(core::magic_name("cl__pathname_device"),true,&core::cl__pathname_device,"(pname &key ((:case scase) :local))");
}
NOINLINE void expose_function_434_helper() {
  expose_function(core::magic_name("cl__pathname_host"),true,&core::cl__pathname_host,"(pname &key ((:case scase) :local))");
}
NOINLINE void expose_function_435_helper() {
  expose_function(core::magic_name("cl__make_pathname"),true,&core::cl__make_pathname,"(&key (host nil hostp) (device nil devicep) (directory nil directoryp) (name nil namep) (type nil typep) (version nil versionp) ((:case scase) :local) (defaults (core:safe-default-pathname-defaults-host-only)))");
}
NOINLINE void expose_function_436_helper() {
  expose_function(core::magic_name("cl__parse_namestring"),true,&core::cl__parse_namestring,"(thing &optional host (defaults (core:safe-default-pathname-defaults)) &key (start 0) end junk-allowed)");
}
NOINLINE void expose_function_437_helper() {
  expose_function(core::magic_name("cl__namestring"),true,&core::cl__namestring,"(pathname)");
}
NOINLINE void expose_function_438_helper() {
  expose_function(core::magic_name("core__coerce_to_filename"),true,&core::core__coerce_to_filename,"(pathname-orig)");
}
NOINLINE void expose_function_439_helper() {
  expose_function(core::magic_name("core__coerce_to_physical_pathname"),true,&core::core__coerce_to_physical_pathname,"(x)");
}
NOINLINE void expose_function_440_helper() {
  expose_function(core::magic_name("core__coerce_to_file_pathname"),true,&core::core__coerce_to_file_pathname,"(tpathname)");
}
NOINLINE void expose_function_441_helper() {
  expose_function(core::magic_name("cl__wild_pathname_p"),true,&core::cl__wild_pathname_p,"(pathname &optional component)");
}
NOINLINE void expose_function_442_helper() {
  expose_function(core::magic_name("cl__merge_pathnames"),true,&core::cl__merge_pathnames,"(arg &optional (default-pathname (core:safe-default-pathname-defaults)) (default-version :newest))");
}
NOINLINE void expose_function_443_helper() {
  expose_function(core::magic_name("cl__pathname"),true,&core::cl__pathname,"");
}
NOINLINE void expose_function_444_helper() {
  expose_function(core::magic_name("core__safe_default_pathname_defaults_host_only"),true,&core::core__safe_default_pathname_defaults_host_only,"");
}
NOINLINE void expose_function_445_helper() {
  expose_function(core::magic_name("core__safe_default_pathname_defaults"),true,&core::core__safe_default_pathname_defaults,"");
}
NOINLINE void expose_function_446_helper() {
  expose_function(core::magic_name("core__float_to_digits"),true,&core::core__float_to_digits,"(digits number position relativep)");
}
NOINLINE void expose_function_447_helper() {
  expose_function(core::magic_name("core__sharp_vertical_bar"),true,&core::core__sharp_vertical_bar,"(stream ch num)");
}
NOINLINE void expose_function_448_helper() {
  expose_function(core::magic_name("core__sharp_minus"),true,&core::core__sharp_minus,"(stream ch num)");
}
NOINLINE void expose_function_449_helper() {
  expose_function(core::magic_name("core__sharp_plus"),true,&core::core__sharp_plus,"(stream ch num)");
}
NOINLINE void expose_function_450_helper() {
  expose_function(core::magic_name("core__reader_feature_p"),true,&core::core__reader_feature_p,"(feature-test)");
}
NOINLINE void expose_function_451_helper() {
  expose_function(core::magic_name("core__sharp_p"),true,&core::core__sharp_p,"(stream ch num)");
}
NOINLINE void expose_function_452_helper() {
  expose_function(core::magic_name("core__sharp_s"),true,&core::core__sharp_s,"(stream ch num)");
}
NOINLINE void expose_function_453_helper() {
  expose_function(core::magic_name("core__sharp_a"),true,&core::core__sharp_a,"(stream ch num)");
}
NOINLINE void expose_function_454_helper() {
  expose_function(core::magic_name("core__sharp_c"),true,&core::core__sharp_c,"(stream ch num)");
}
NOINLINE void expose_function_455_helper() {
  expose_function(core::magic_name("core__sharp_x"),true,&core::core__sharp_x,"(stream ch num)");
}
NOINLINE void expose_function_456_helper() {
  expose_function(core::magic_name("core__sharp_o"),true,&core::core__sharp_o,"(stream ch num)");
}
NOINLINE void expose_function_457_helper() {
  expose_function(core::magic_name("core__sharp_b"),true,&core::core__sharp_b,"(stream ch num)");
}
NOINLINE void expose_function_458_helper() {
  expose_function(core::magic_name("core__sharp_r"),true,&core::core__sharp_r,"(stream subchar radix)");
}
NOINLINE void expose_function_459_helper() {
  expose_function(core::magic_name("core__sharp_colon"),true,&core::core__sharp_colon,"(stream ch num)");
}
NOINLINE void expose_function_460_helper() {
  expose_function(core::magic_name("core__sharp_asterisk"),true,&core::core__sharp_asterisk,"(stream ch num)");
}
NOINLINE void expose_function_461_helper() {
  expose_function(core::magic_name("core__sharp_left_parenthesis"),true,&core::core__sharp_left_parenthesis,"(stream ch num)");
}
NOINLINE void expose_function_462_helper() {
  expose_function(core::magic_name("core__sharp_single_quote"),true,&core::core__sharp_single_quote,"(stream ch num)");
}
NOINLINE void expose_function_463_helper() {
  expose_function(core::magic_name("core__sharp_dot"),true,&core::core__sharp_dot,"(stream ch num)");
}
NOINLINE void expose_function_464_helper() {
  expose_function(core::magic_name("core__sharp_backslash"),true,&core::core__sharp_backslash,"(stream ch num)");
}
NOINLINE void expose_function_465_helper() {
  expose_function(core::magic_name("core__dispatch_macro_character"),true,&core::core__dispatch_macro_character,"(sin ch)");
}
NOINLINE void expose_function_466_helper() {
  expose_function(core::magic_name("core__reader_skip_semicolon_comment"),true,&core::core__reader_skip_semicolon_comment,"(sin ch)");
}
NOINLINE void expose_function_467_helper() {
  expose_function(core::magic_name("core__reader_quote"),true,&core::core__reader_quote,"(sin ch)");
}
NOINLINE void expose_function_468_helper() {
  expose_function(core::magic_name("core__reader_error_unmatched_close_parenthesis"),true,&core::core__reader_error_unmatched_close_parenthesis,"(sin ch)");
}
NOINLINE void expose_function_469_helper() {
  expose_function(core::magic_name("core__reader_list_allow_consing_dot"),true,&core::core__reader_list_allow_consing_dot,"(sin ch)");
}
NOINLINE void expose_function_470_helper() {
  expose_function(core::magic_name("core__reader_comma_form"),true,&core::core__reader_comma_form,"(sin ch)");
}
NOINLINE void expose_function_471_helper() {
  expose_function(core::magic_name("core__reader_backquoted_expression"),true,&core::core__reader_backquoted_expression,"(sin ch)");
}
NOINLINE void expose_function_472_helper() {
  expose_function(core::magic_name("core__reader_double_quote_string"),true,&core::core__reader_double_quote_string,"(stream chr)");
}
NOINLINE void expose_function_473_helper() {
  expose_function(core::magic_name("cl__set_macro_character"),true,&core::cl__set_macro_character,"(ch func_desig &optional non-terminating-p (readtable *readtable*))");
}
NOINLINE void expose_function_474_helper() {
  expose_function(core::magic_name("cl__get_dispatch_macro_character"),true,&core::cl__get_dispatch_macro_character,"(dispChar subChar &optional (readtable *readtable*))");
}
NOINLINE void expose_function_475_helper() {
  expose_function(core::magic_name("cl__set_dispatch_macro_character"),true,&core::cl__set_dispatch_macro_character,"(dispChar subChar newFunction &optional (readtable *readtable*))");
}
NOINLINE void expose_function_476_helper() {
  expose_function(core::magic_name("core__readtable_case_set"),true,&core::core__readtable_case_set,"(readtable mode)");
}
NOINLINE void expose_function_477_helper() {
  expose_function(core::magic_name("cl__readtable_case"),true,&core::cl__readtable_case,"(readtable)");
}
NOINLINE void expose_function_478_helper() {
  expose_function(core::magic_name("cl__copy_readtable"),true,&core::cl__copy_readtable,"(&optional (from-readtable cl:*readtable*) to-readtable)");
}
NOINLINE void expose_function_479_helper() {
  expose_function(core::magic_name("cl__get_macro_character"),true,&core::cl__get_macro_character,"(char &optional readtable)");
}
NOINLINE void expose_function_480_helper() {
  expose_function(core::magic_name("cl__make_dispatch_macro_character"),true,&core::cl__make_dispatch_macro_character,"(char &optional non-terminating-p (readtable *readtable*))");
}
NOINLINE void expose_function_481_helper() {
  expose_function(core::magic_name("cl__set_syntax_from_char"),true,&core::cl__set_syntax_from_char,"(tochar fromchar &optional (toreadtable *readtable*) (fromreadtable nil fromreadtablep))");
}
NOINLINE void expose_function_482_helper() {
  expose_function(core::magic_name("core__make_branch_snode"),true,&core::core__make_branch_snode,"");
}
NOINLINE void expose_function_483_helper() {
  expose_function(core::magic_name("core__pointer_address"),true,&core::core__pointer_address,"(arg)");
}
NOINLINE void expose_function_484_helper() {
  expose_function(core::magic_name("core__pointer_delete"),true,&core::core__pointer_delete,"(arg)");
}
NOINLINE void expose_function_485_helper() {
  expose_function(core::magic_name("core__pointer_release"),true,&core::core__pointer_release,"(arg)");
}
NOINLINE void expose_function_486_helper() {
  expose_function(core::magic_name("core__make_regex"),true,&core::core__make_regex,"(regex-str)");
}
NOINLINE void expose_function_487_helper() {
  expose_function(core::magic_name("core__write_ugly_object"),true,&core::core__write_ugly_object,"(obj &optional strm)");
}
NOINLINE void expose_function_488_helper() {
  expose_function(core::magic_name("core__test_profile_timer"),true,&core::core__test_profile_timer,"(delay)");
}
NOINLINE void expose_function_489_helper() {
  expose_function(core::magic_name("core__clock_gettime_nanoseconds"),true,&core::core__clock_gettime_nanoseconds,"");
}
NOINLINE void expose_function_490_helper() {
  expose_function(core::magic_name("core__deep_copy"),true,&core::core__deep_copy,"(obj)");
}
NOINLINE void expose_function_491_helper() {
  expose_function(core::magic_name("core__sl_boundp"),true,&core::core__sl_boundp,"(arg)");
}
NOINLINE void expose_function_492_helper() {
  expose_function(core::magic_name("core__decode"),true,&core::core__decode,"(obj arg)");
}
NOINLINE void expose_function_493_helper() {
  expose_function(core::magic_name("core__encode"),true,&core::core__encode,"(arg)");
}
NOINLINE void expose_function_494_helper() {
  expose_function(core::magic_name("core__is_nil"),true,&core::core__is_nil,"(arg)");
}
NOINLINE void expose_function_495_helper() {
  expose_function(core::magic_name("core__instancep"),true,&core::core__instancep,"(obj)");
}
NOINLINE void expose_function_496_helper() {
  expose_function(core::magic_name("core__instance_ref"),true,&core::core__instance_ref,"(obj idx)");
}
NOINLINE void expose_function_497_helper() {
  expose_function(core::magic_name("core__instance_set"),true,&core::core__instance_set,"(obj idx val)");
}
NOINLINE void expose_function_498_helper() {
  expose_function(core::magic_name("core__instance_sig_set"),true,&core::core__instance_sig_set,"(arg)");
}
NOINLINE void expose_function_499_helper() {
  expose_function(core::magic_name("core__instance_sig"),true,&core::core__instance_sig,"(arg)");
}
NOINLINE void expose_function_500_helper() {
  expose_function(core::magic_name("core__class_name_as_string"),true,&core::core__class_name_as_string,"(arg)");
}
NOINLINE void expose_function_501_helper() {
  expose_function(core::magic_name("core__instance_class"),true,&core::core__instance_class,"(arg)");
}
NOINLINE void expose_function_502_helper() {
  expose_function(core::magic_name("core__implementation_class"),true,&core::core__implementation_class,"(arg)");
}
NOINLINE void expose_function_503_helper() {
  expose_function(core::magic_name("cl__copy_tree"),true,&core::cl__copy_tree,"(arg)");
}
NOINLINE void expose_function_504_helper() {
  expose_function(core::magic_name("core__low_level_describe"),true,&core::core__low_level_describe,"(arg)");
}
NOINLINE void expose_function_505_helper() {
  expose_function(core::magic_name("core__print_cxx_object"),true,&core::core__print_cxx_object,"(obj stream)");
}
NOINLINE void expose_function_506_helper() {
  expose_function(core::magic_name("core__fieldsp"),true,&core::core__fieldsp,"(obj)");
}
NOINLINE void expose_function_507_helper() {
  expose_function(core::magic_name("core__make_cxx_object"),true,&core::core__make_cxx_object,"(class-name &rest args)");
}
NOINLINE void expose_function_508_helper() {
  expose_function(core::magic_name("core__all_source_files"),true,&core::core__all_source_files,"");
}
NOINLINE void expose_function_509_helper() {
  expose_function(core::magic_name("CORE:forget_all_single_dispatch_generic_functions"),true,&core::Lisp_O::forget_all_single_dispatch_generic_functions,"");
}
NOINLINE void expose_function_510_helper() {
  expose_function(core::magic_name("CORE:setf_find_single_dispatch_generic_function"),true,&core::Lisp_O::setf_find_single_dispatch_generic_function,"(gf-symbol gf)");
}
NOINLINE void expose_function_511_helper() {
  expose_function(core::magic_name("CORE:find_single_dispatch_generic_function"),true,&core::Lisp_O::find_single_dispatch_generic_function,"(gf-symbol &optional errorp)");
}
NOINLINE void expose_function_512_helper() {
  expose_function(core::magic_name("cl__not"),true,&core::cl__not,"(arg)");
}
NOINLINE void expose_function_513_helper() {
  expose_function(core::magic_name("core__repr"),true,&core::core__repr,"(arg)");
}
NOINLINE void expose_function_514_helper() {
  expose_function(core::magic_name("core__is_sub_class_of"),true,&core::core__is_sub_class_of,"(tag mc)");
}
NOINLINE void expose_function_515_helper() {
  expose_function(core::magic_name("core__is_assignable_to"),true,&core::core__is_assignable_to,"(tag secondArgument)");
}
NOINLINE void expose_function_516_helper() {
  expose_function(core::magic_name("cl__cerror"),true,&core::cl__cerror,"(cformat eformat &rest arguments)");
}
NOINLINE void expose_function_517_helper() {
  expose_function(core::magic_name("cl__error"),true,&core::cl__error,"(datum &rest arguments)");
}
NOINLINE void expose_function_518_helper() {
  expose_function(core::magic_name("core__invoke_internal_debugger_from_gdb"),true,&core::core__invoke_internal_debugger_from_gdb,"");
}
NOINLINE void expose_function_519_helper() {
  expose_function(core::magic_name("core__single_dispatch_generic_function_table"),true,&core::core__single_dispatch_generic_function_table,"");
}
NOINLINE void expose_function_520_helper() {
  expose_function(core::magic_name("core__invoke_internal_debugger"),true,&core::core__invoke_internal_debugger,"(&optional condition)");
}
NOINLINE void expose_function_521_helper() {
  expose_function(core::magic_name("core__universal_error_handler"),true,&core::core__universal_error_handler,"(continue-string datum initializers)");
}
NOINLINE void expose_function_522_helper() {
  expose_function(core::magic_name("cl__intern"),true,&core::cl__intern,"(symbol_name &optional (package-desig *package*))");
}
NOINLINE void expose_function_523_helper() {
  expose_function(core::magic_name("core__export_to_python"),true,&core::core__export_to_python,"(symbolsDesig)");
}
NOINLINE void expose_function_524_helper() {
  expose_function(core::magic_name("cl__export"),true,&core::cl__export,"(symDes &optional (packageDes *package*))");
}
NOINLINE void expose_function_525_helper() {
  expose_function(core::magic_name("core__debug_log_off"),true,&core::core__debug_log_off,"");
}
NOINLINE void expose_function_526_helper() {
  expose_function(core::magic_name("core__debug_log_on"),true,&core::core__debug_log_on,"");
}
NOINLINE void expose_function_527_helper() {
  expose_function(core::magic_name("core__is_top_level_script"),true,&core::core__is_top_level_script,"");
}
NOINLINE void expose_function_528_helper() {
  expose_function(core::magic_name("core__set_current_working_directory"),true,&core::core__set_current_working_directory,"(dir)");
}
NOINLINE void expose_function_529_helper() {
  expose_function(core::magic_name("core__database_dir"),true,&core::core__database_dir,"");
}
NOINLINE void expose_function_530_helper() {
  expose_function(core::magic_name("core__lisp_code_path"),true,&core::core__lisp_code_path,"(&optional rel-path)");
}
NOINLINE void expose_function_531_helper() {
  expose_function(core::magic_name("core__library_path"),true,&core::core__library_path,"(&optional rel-path)");
}
NOINLINE void expose_function_532_helper() {
  expose_function(core::magic_name("core__script_dir"),true,&core::core__script_dir,"");
}
NOINLINE void expose_function_533_helper() {
  expose_function(core::magic_name("core__source_line_column"),true,&core::core__source_line_column,"");
}
NOINLINE void expose_function_534_helper() {
  expose_function(core::magic_name("core__source_file_name"),true,&core::core__source_file_name,"");
}
NOINLINE void expose_function_535_helper() {
  expose_function(core::magic_name("cl__sort"),true,&core::cl__sort,"(sequence predicate)");
}
NOINLINE void expose_function_536_helper() {
  expose_function(core::magic_name("core__sorted"),true,&core::core__sorted,"(unsorted)");
}
NOINLINE void expose_function_537_helper() {
  expose_function(core::magic_name("cl__apropos"),true,&core::cl__apropos,"(string_desig &optional package_desig)");
}
NOINLINE void expose_function_538_helper() {
  expose_function(core::magic_name("cl__macroexpand"),true,&core::cl__macroexpand,"(form &optional env)");
}
NOINLINE void expose_function_539_helper() {
  expose_function(core::magic_name("cl__macroexpand_1"),true,&core::cl__macroexpand_1,"(form &optional env)");
}
NOINLINE void expose_function_540_helper() {
  expose_function(core::magic_name("core__mpi_size"),true,&core::core__mpi_size,"");
}
NOINLINE void expose_function_541_helper() {
  expose_function(core::magic_name("core__mpi_rank"),true,&core::core__mpi_rank,"");
}
NOINLINE void expose_function_542_helper() {
  expose_function(core::magic_name("core__mpi_enabled"),true,&core::core__mpi_enabled,"");
}
NOINLINE void expose_function_543_helper() {
  expose_function(core::magic_name("core__select_package"),true,&core::core__select_package,"(package-designator)");
}
NOINLINE void expose_function_544_helper() {
  expose_function(core::magic_name("cl__find_package"),true,&core::cl__find_package,"(name_desig)");
}
NOINLINE void expose_function_545_helper() {
  expose_function(core::magic_name("core__find_file_in_lisp_path"),true,&core::core__find_file_in_lisp_path,"(partialPath)");
}
NOINLINE void expose_function_546_helper() {
  expose_function(core::magic_name("core__setf_find_class"),true,&core::core__setf_find_class,"(new-value name)");
}
NOINLINE void expose_function_547_helper() {
  expose_function(core::magic_name("cl__find_class"),true,&core::cl__find_class,"(symbol &optional (errorp t) environment)");
}
NOINLINE void expose_function_548_helper() {
  expose_function(core::magic_name("core__getline"),true,&core::core__getline,"(&optional (prompt \"\"))");
}
NOINLINE void expose_function_549_helper() {
  expose_function(core::magic_name("core__member1"),true,&core::core__member1,"(item list test test-not key)");
}
NOINLINE void expose_function_550_helper() {
  expose_function(core::magic_name("core__member_test"),true,&core::core__member_test,"(item list &key key test test-not)");
}
NOINLINE void expose_function_551_helper() {
  expose_function(core::magic_name("cl__member"),true,&core::cl__member,"(item list &key key test test-not)");
}
NOINLINE void expose_function_552_helper() {
  expose_function(core::magic_name("cl__assoc"),true,&core::cl__assoc,"(item alist &key key test test-not)");
}
NOINLINE void expose_function_553_helper() {
  expose_function(core::magic_name("cl__acons"),true,&core::cl__acons,"(key datum alist)");
}
NOINLINE void expose_function_554_helper() {
  expose_function(core::magic_name("core__quit"),true,&core::core__quit,"(&optional (exit-value 0))");
}
NOINLINE void expose_function_555_helper() {
  expose_function(core::magic_name("core__exit"),true,&core::core__exit,"(&optional (exit-value 0))");
}
NOINLINE void expose_function_556_helper() {
  expose_function(core::magic_name("core__setup_stack_monitor"),true,&core::core__setup_stack_monitor,"(&key warn-size sample-size)");
}
NOINLINE void expose_function_557_helper() {
  expose_function(core::magic_name("core__stack_limit"),true,&core::core__stack_limit,"");
}
NOINLINE void expose_function_558_helper() {
  expose_function(core::magic_name("core__stack_monitor"),true,&core::core__stack_monitor,"");
}
NOINLINE void expose_function_559_helper() {
  expose_function(core::magic_name("core__stack_used"),true,&core::core__stack_used,"");
}
NOINLINE void expose_function_560_helper() {
  expose_function(core::magic_name("core__low_level_repl"),true,&core::core__low_level_repl,"");
}
NOINLINE void expose_function_561_helper() {
  expose_function(core::magic_name("core__list_of_all_special_operators"),true,&core::core__list_of_all_special_operators,"");
}
NOINLINE void expose_function_562_helper() {
  expose_function(core::magic_name("cl__boole"),true,&core::cl__boole,"(op arg1 arg2)");
}
NOINLINE void expose_function_563_helper() {
  expose_function(core::magic_name("cl__logbitp"),true,&core::cl__logbitp,"");
}
NOINLINE void expose_function_564_helper() {
  expose_function(core::magic_name("core__bit_array_op"),true,&core::core__bit_array_op,"(op x y &optional r)");
}
NOINLINE void expose_function_565_helper() {
  expose_function(core::magic_name("cl__write_sequence"),true,&core::cl__write_sequence,"(seq stream &key (start 0) end)");
}
NOINLINE void expose_function_566_helper() {
  expose_function(core::magic_name("core__file_column"),true,&core::core__file_column,"(arg)");
}
NOINLINE void expose_function_567_helper() {
  expose_function(core::magic_name("cl__unread_char"),true,&core::cl__unread_char,"(char &optional strm)");
}
NOINLINE void expose_function_568_helper() {
  expose_function(core::magic_name("cl__finish_output"),true,&core::cl__finish_output,"(&optional strm)");
}
NOINLINE void expose_function_569_helper() {
  expose_function(core::magic_name("cl__force_output"),true,&core::cl__force_output,"(&optional strm)");
}
NOINLINE void expose_function_570_helper() {
  expose_function(core::magic_name("cl__listen"),true,&core::cl__listen,"(&optional dstrm)");
}
NOINLINE void expose_function_571_helper() {
  expose_function(core::magic_name("cl__clear_output"),true,&core::cl__clear_output,"(&optional dstrm)");
}
NOINLINE void expose_function_572_helper() {
  expose_function(core::magic_name("cl__clear_input"),true,&core::cl__clear_input,"(&optional dstrm)");
}
NOINLINE void expose_function_573_helper() {
  expose_function(core::magic_name("cl__write_char"),true,&core::cl__write_char,"(string &optional output-stream)");
}
NOINLINE void expose_function_574_helper() {
  expose_function(core::magic_name("cl__write_byte"),true,&core::cl__write_byte,"(byte &optional output-stream)");
}
NOINLINE void expose_function_575_helper() {
  expose_function(core::magic_name("cl__write_line"),true,&core::cl__write_line,"(string &optional output-stream &key (start 0) end)");
}
NOINLINE void expose_function_576_helper() {
  expose_function(core::magic_name("cl__write_string"),true,&core::cl__write_string,"(string &optional output-stream &key (start 0) end)");
}
NOINLINE void expose_function_577_helper() {
  expose_function(core::magic_name("cl__fresh_line"),true,&core::cl__fresh_line,"(&optional outputStream)");
}
NOINLINE void expose_function_578_helper() {
  expose_function(core::magic_name("cl__terpri"),true,&core::cl__terpri,"(&optional output-stream)");
}
NOINLINE void expose_function_579_helper() {
  expose_function(core::magic_name("cl__read_line"),true,&core::cl__read_line,"(&optional input-stream (eof-error-p t) eof-value recursive-p)");
}
NOINLINE void expose_function_580_helper() {
  expose_function(core::magic_name("cl__read_from_string"),true,&core::cl__read_from_string,"(content &optional (eof-error-p t) eof-value &key (start 0) end preserve-whitespace)");
}
NOINLINE void expose_function_581_helper() {
  expose_function(core::magic_name("cl__read_char_no_hang"),true,&core::cl__read_char_no_hang,"(&optional strm (eof_error_p t) eof_value recursive_p)");
}
NOINLINE void expose_function_582_helper() {
  expose_function(core::magic_name("cl__read_char"),true,&core::cl__read_char,"(&optional strm (eof_error_p t) eof_value recursive_p)");
}
NOINLINE void expose_function_583_helper() {
  expose_function(core::magic_name("cl__peek_char"),true,&core::cl__peek_char,"(&optional peek_type strm (eof_errorp t) eof_value recursivep)");
}
NOINLINE void expose_function_584_helper() {
  expose_function(core::magic_name("cl__read_byte"),true,&core::cl__read_byte,"(&optional strm (eof_error_p t) eof_value)");
}
NOINLINE void expose_function_585_helper() {
  expose_function(core::magic_name("core__input_stream_source_pos_info"),true,&core::core__input_stream_source_pos_info,"(arg)");
}
NOINLINE void expose_function_586_helper() {
  expose_function(core::magic_name("core__stream_column"),true,&core::core__stream_column,"(stream)");
}
NOINLINE void expose_function_587_helper() {
  expose_function(core::magic_name("core__stream_linenumber"),true,&core::core__stream_linenumber,"(stream)");
}
NOINLINE void expose_function_588_helper() {
  expose_function(core::magic_name("cl__close"),true,&core::cl__close,"(strm &key abort)");
}
NOINLINE void expose_function_589_helper() {
  expose_function(core::magic_name("cl__open"),true,&core::cl__open,"(filename &key (direction :input) (element-type 'base-char) (if-exists nil iesp) (if-does-not-exist nil idnesp) (external-format :default) (cstream T))");
}
NOINLINE void expose_function_590_helper() {
  expose_function(core::magic_name("cl__streamp"),true,&core::cl__streamp,"(arg)");
}
NOINLINE void expose_function_591_helper() {
  expose_function(core::magic_name("cl__interactive_stream_p"),true,&core::cl__interactive_stream_p,"(arg)");
}
NOINLINE void expose_function_592_helper() {
  expose_function(core::magic_name("cl__output_stream_p"),true,&core::cl__output_stream_p,"(arg)");
}
NOINLINE void expose_function_593_helper() {
  expose_function(core::magic_name("cl__input_stream_p"),true,&core::cl__input_stream_p,"(strm)");
}
NOINLINE void expose_function_594_helper() {
  expose_function(core::magic_name("cl__file_position"),true,&core::cl__file_position,"(file-stream &optional position)");
}
NOINLINE void expose_function_595_helper() {
  expose_function(core::magic_name("cl__file_length"),true,&core::cl__file_length,"(strm)");
}
NOINLINE void expose_function_596_helper() {
  expose_function(core::magic_name("cl__read_sequence"),true,&core::cl__read_sequence,"(sequence stream &key (start 0) end)");
}
NOINLINE void expose_function_597_helper() {
  expose_function(core::magic_name("core__do_write_sequence"),true,&core::core__do_write_sequence,"(seq stream start end)");
}
NOINLINE void expose_function_598_helper() {
  expose_function(core::magic_name("core__file_stream_fd"),true,&core::core__file_stream_fd,"(s)");
}
NOINLINE void expose_function_599_helper() {
  expose_function(core::magic_name("core__set_buffering_mode"),true,&core::core__set_buffering_mode,"(stream mode)");
}
NOINLINE void expose_function_600_helper() {
  expose_function(core::magic_name("cl__synonym_stream_symbol"),true,&core::cl__synonym_stream_symbol,"(s)");
}
NOINLINE void expose_function_601_helper() {
  expose_function(core::magic_name("cl__make_synonym_stream"),true,&core::cl__make_synonym_stream,"(strm1)");
}
NOINLINE void expose_function_602_helper() {
  expose_function(core::magic_name("cl__make_concatenated_stream"),true,&core::cl__make_concatenated_stream,"(&rest ap)");
}
NOINLINE void expose_function_603_helper() {
  expose_function(core::magic_name("cl__make_echo_stream"),true,&core::cl__make_echo_stream,"(strm1 strm2)");
}
NOINLINE void expose_function_604_helper() {
  expose_function(core::magic_name("cl__make_broadcast_stream"),true,&core::cl__make_broadcast_stream,"(&rest ap)");
}
NOINLINE void expose_function_605_helper() {
  expose_function(core::magic_name("cl__two_way_stream_output_stream"),true,&core::cl__two_way_stream_output_stream,"(strm)");
}
NOINLINE void expose_function_606_helper() {
  expose_function(core::magic_name("cl__two_way_stream_input_stream"),true,&core::cl__two_way_stream_input_stream,"(strm)");
}
NOINLINE void expose_function_607_helper() {
  expose_function(core::magic_name("cl__make_two_way_stream"),true,&core::cl__make_two_way_stream,"(istrm ostrm)");
}
NOINLINE void expose_function_608_helper() {
  expose_function(core::magic_name("cl__make_string_input_stream"),true,&core::cl__make_string_input_stream,"(strng &optional (istart 0) iend)");
}
NOINLINE void expose_function_609_helper() {
  expose_function(core::magic_name("cl__get_output_stream_string"),true,&core::cl__get_output_stream_string,"(strm)");
}
NOINLINE void expose_function_610_helper() {
  expose_function(core::magic_name("cl__make_string_output_stream"),true,&core::cl__make_string_output_stream,"(&key (element-type 'character))");
}
NOINLINE void expose_function_611_helper() {
  expose_function(core::magic_name("core__make_string_output_stream_from_string"),true,&core::core__make_string_output_stream_from_string,"(s)");
}
NOINLINE void expose_function_612_helper() {
  expose_function(core::magic_name("CORE:makeLambdaListHandler"),true,&core::LambdaListHandler_O::makeLambdaListHandler,"(lambda-list &optional declares (context 'core::function))");
}
NOINLINE void expose_function_613_helper() {
  expose_function(core::magic_name("core__process_lambda_list"),true,&core::core__process_lambda_list,"(vl context)");
}
NOINLINE void expose_function_614_helper() {
  expose_function(core::magic_name("core__process_single_dispatch_lambda_list"),true,&core::core__process_single_dispatch_lambda_list,"(lambda-list)");
}
NOINLINE void expose_function_615_helper() {
  expose_function(core::magic_name("core__process_macro_lambda_list"),true,&core::core__process_macro_lambda_list,"(lambda-list)");
}
NOINLINE void expose_function_616_helper() {
  expose_function(core::magic_name("core__magic_name"),true,&core::core__magic_name,"(name)");
}
NOINLINE void expose_function_617_helper() {
  expose_function(core::magic_name("core__lispify_name"),true,&core::core__lispify_name,"(name)");
}
NOINLINE void expose_function_618_helper() {
  expose_function(core::magic_name("CORE:makeStackValueEnvironment"),true,&core::StackValueEnvironment_O::make,"");
}
NOINLINE void expose_function_619_helper() {
  expose_function(core::magic_name("CORE:makeSymbolMacroletEnvironment"),true,&core::SymbolMacroletEnvironment_O::make,"");
}
NOINLINE void expose_function_620_helper() {
  expose_function(core::magic_name("CORE:makeMacroletEnvironment"),true,&core::MacroletEnvironment_O::make,"");
}
NOINLINE void expose_function_621_helper() {
  expose_function(core::magic_name("CORE:makeTagbodyEnvironment"),true,&core::TagbodyEnvironment_O::make,"");
}
NOINLINE void expose_function_622_helper() {
  expose_function(core::magic_name("CORE:makeFunctionContainerEnvironment"),true,&core::FunctionContainerEnvironment_O::make,"");
}
NOINLINE void expose_function_623_helper() {
  expose_function(core::magic_name("CORE:makeCatchEnvironment"),true,&core::CatchEnvironment_O::make,"");
}
NOINLINE void expose_function_624_helper() {
  expose_function(core::magic_name("CORE:makeBlockEnvironment"),true,&core::BlockEnvironment_O::make,"");
}
NOINLINE void expose_function_625_helper() {
  expose_function(core::magic_name("CORE:makeUnwindProtectEnvironment"),true,&core::UnwindProtectEnvironment_O::make,"");
}
NOINLINE void expose_function_626_helper() {
  expose_function(core::magic_name("CORE:makeFunctionValueEnvironment"),true,&core::FunctionValueEnvironment_O::createForEntries,"");
}
NOINLINE void expose_function_627_helper() {
  expose_function(core::magic_name("CORE:makeValueEnvironmentForLocallySpecialEntries"),true,&core::ValueEnvironment_O::createForLocallySpecialEntries,"");
}
NOINLINE void expose_function_628_helper() {
  expose_function(core::magic_name("CORE:makeValueEnvironmentForNumberOfEntries"),true,&core::ValueEnvironment_O::createForNumberOfEntries,"");
}
NOINLINE void expose_function_629_helper() {
  expose_function(core::magic_name("CORE:makeValueEnvironment"),true,&core::ValueEnvironment_O::createForLambdaListHandler,"");
}
NOINLINE void expose_function_630_helper() {
  expose_function("CORE:CURRENT-VISIBLE-ENVIRONMENT",true,&core::Environment_O::clasp_currentVisibleEnvironment,"");
}
NOINLINE void expose_function_631_helper() {
  expose_function(core::magic_name("core__environment_id"),true,&core::core__environment_id,"(env)");
}
NOINLINE void expose_function_632_helper() {
  expose_function(core::magic_name("core__runtime_environment"),true,&core::core__runtime_environment,"(env)");
}
NOINLINE void expose_function_633_helper() {
  expose_function(core::magic_name("core__environment_type_list"),true,&core::core__environment_type_list,"(env)");
}
NOINLINE void expose_function_634_helper() {
  expose_function(core::magic_name("core__environment_list"),true,&core::core__environment_list,"(env)");
}
NOINLINE void expose_function_635_helper() {
  expose_function(core::magic_name("core__environment_activation_frame"),true,&core::core__environment_activation_frame,"(env)");
}
NOINLINE void expose_function_636_helper() {
  expose_function(core::magic_name("core__lexical_macro_function"),true,&core::core__lexical_macro_function,"(name env)");
}
NOINLINE void expose_function_637_helper() {
  expose_function(core::magic_name("core__lexical_function"),true,&core::core__lexical_function,"(name env)");
}
NOINLINE void expose_function_638_helper() {
  expose_function(core::magic_name("core__environment_debug_values"),true,&core::core__environment_debug_values,"(frame)");
}
NOINLINE void expose_function_639_helper() {
  expose_function(core::magic_name("core__environment_debug_names"),true,&core::core__environment_debug_names,"(frame)");
}
NOINLINE void expose_function_640_helper() {
  expose_function(core::magic_name("core__environment_length"),true,&core::core__environment_length,"(frame)");
}
NOINLINE void expose_function_641_helper() {
  expose_function(core::magic_name("core__classify_return_from_symbol"),true,&core::core__classify_return_from_symbol,"(env sym)");
}
NOINLINE void expose_function_642_helper() {
  expose_function(core::magic_name("cl__nreconc"),true,&core::cl__nreconc,"(list tail)");
}
NOINLINE void expose_function_643_helper() {
  expose_function(core::magic_name("cl__revappend"),true,&core::cl__revappend,"(list tail)");
}
NOINLINE void expose_function_644_helper() {
  expose_function(core::magic_name("cl__nconc"),true,&core::cl__nconc,"(&rest lists)");
}
NOINLINE void expose_function_645_helper() {
  expose_function(core::magic_name("cl__last"),true,&core::cl__last,"(list &optional (on 1))");
}
NOINLINE void expose_function_646_helper() {
  expose_function(core::magic_name("cl__listSTAR"),true,&core::cl__listSTAR,"(&rest objects)");
}
NOINLINE void expose_function_647_helper() {
  expose_function(core::magic_name("cl__list"),true,&core::cl__list,"(&rest objects)");
}
NOINLINE void expose_function_648_helper() {
  expose_function(core::magic_name("cl__nbutlast"),true,&core::cl__nbutlast,"(list &optional (n 1))");
}
NOINLINE void expose_function_649_helper() {
  expose_function(core::magic_name("cl__butlast"),true,&core::cl__butlast,"(list &optional (n 1))");
}
NOINLINE void expose_function_650_helper() {
  expose_function(core::magic_name("cl__copy_list"),true,&core::cl__copy_list,"(arg)");
}
NOINLINE void expose_function_651_helper() {
  expose_function(core::magic_name("cl__nthcdr"),true,&core::cl__nthcdr,"(idx arg)");
}
NOINLINE void expose_function_652_helper() {
  expose_function(core::magic_name("cl__nth"),true,&core::cl__nth,"(idx arg)");
}
NOINLINE void expose_function_653_helper() {
  expose_function(core::magic_name("cl__rassoc"),true,&core::cl__rassoc,"(item a-list &key test test-not key)");
}
NOINLINE void expose_function_654_helper() {
  expose_function(core::magic_name("core__function_source_code"),true,&core::core__function_source_code,"(fn)");
}
NOINLINE void expose_function_655_helper() {
  expose_function(core::magic_name("cl__function_lambda_expression"),true,&core::cl__function_lambda_expression,"(fn)");
}
NOINLINE void expose_function_656_helper() {
  expose_function(core::magic_name("core__set_kind"),true,&core::core__set_kind,"(fn kind)");
}
NOINLINE void expose_function_657_helper() {
  expose_function(core::magic_name("core__function_source_pos_info"),true,&core::core__function_source_pos_info,"(function)");
}
NOINLINE void expose_function_658_helper() {
  expose_function(core::magic_name("core__function_lambda_list"),true,&core::core__function_lambda_list,"(function)");
}
NOINLINE void expose_function_659_helper() {
  expose_function(core::magic_name("core__mkdir"),true,&core::core__mkdir,"(dir mode)");
}
NOINLINE void expose_function_660_helper() {
  expose_function(core::magic_name("core__unix_get_local_time_zone"),true,&core::core__unix_get_local_time_zone,"");
}
NOINLINE void expose_function_661_helper() {
  expose_function(core::magic_name("core__unix_daylight_saving_time"),true,&core::core__unix_daylight_saving_time,"(unix-time)");
}
NOINLINE void expose_function_662_helper() {
  expose_function(core::magic_name("cl__directory"),true,&core::cl__directory,"(mask &key (resolve-symlinks t) &allow-other-keys)");
}
NOINLINE void expose_function_663_helper() {
  expose_function(core::magic_name("core__copy_file"),true,&core::core__copy_file,"(orig dest)");
}
NOINLINE void expose_function_664_helper() {
  expose_function(core::magic_name("core__chmod"),true,&core::core__chmod,"(file mode)");
}
NOINLINE void expose_function_665_helper() {
  expose_function(core::magic_name("core__rmdir"),true,&core::core__rmdir,"(directory)");
}
NOINLINE void expose_function_666_helper() {
  expose_function(core::magic_name("core__mkstemp"),true,&core::core__mkstemp,"(template)");
}
NOINLINE void expose_function_667_helper() {
  expose_function(core::magic_name("cl__user_homedir_pathname"),true,&core::cl__user_homedir_pathname,"(&optional host)");
}
NOINLINE void expose_function_668_helper() {
  expose_function(core::magic_name("cl__file_author"),true,&core::cl__file_author,"(file)");
}
NOINLINE void expose_function_669_helper() {
  expose_function(core::magic_name("cl__file_write_date"),true,&core::cl__file_write_date,"(pathspec)");
}
NOINLINE void expose_function_670_helper() {
  expose_function(core::magic_name("cl__probe_file"),true,&core::cl__probe_file,"(filespec)");
}
NOINLINE void expose_function_671_helper() {
  expose_function(core::magic_name("cl__delete_file"),true,&core::cl__delete_file,"(file)");
}
NOINLINE void expose_function_672_helper() {
  expose_function(core::magic_name("cl__rename_file"),true,&core::cl__rename_file,"(oldn newn &key (if-exists :error))");
}
NOINLINE void expose_function_673_helper() {
  expose_function(core::magic_name("cl__truename"),true,&core::cl__truename,"(orig-pathname)");
}
NOINLINE void expose_function_674_helper() {
  expose_function(core::magic_name("core__file_truename"),true,&core::core__file_truename,"(pathname filename follow-links)");
}
NOINLINE void expose_function_675_helper() {
  expose_function(core::magic_name("core__readlink"),true,&core::core__readlink,"(filename)");
}
NOINLINE void expose_function_676_helper() {
  expose_function(core::magic_name("core__file_kind"),true,&core::core__file_kind,"(filename follow-links)");
}
NOINLINE void expose_function_677_helper() {
  expose_function(core::magic_name("core__current_dir"),true,&core::core__current_dir,"");
}
NOINLINE void expose_function_678_helper() {
  expose_function(core::magic_name("ext__chdir"),true,&core::ext__chdir,"(pathname)");
}
NOINLINE void expose_function_679_helper() {
  expose_function(core::magic_name("core__getppid"),true,&core::core__getppid,"");
}
NOINLINE void expose_function_680_helper() {
  expose_function(core::magic_name("core__getpid"),true,&core::core__getpid,"");
}
NOINLINE void expose_function_681_helper() {
  expose_function(core::magic_name("core__waitpid"),true,&core::core__waitpid,"(pid options)");
}
NOINLINE void expose_function_682_helper() {
  expose_function(core::magic_name("core__fork"),true,&core::core__fork,"");
}
NOINLINE void expose_function_683_helper() {
  expose_function(core::magic_name("CORE:allocateForeignObject"),true,&core::ForeignData_O::allocateForeignObject,"");
}
NOINLINE void expose_function_684_helper() {
  expose_function(core::magic_name("core__make_source_manager"),true,&core::core__make_source_manager,"");
}
NOINLINE void expose_function_685_helper() {
  expose_function(core::magic_name("core__dump_source_manager"),true,&core::core__dump_source_manager,"(dumpAll)");
}
NOINLINE void expose_function_686_helper() {
  expose_function(core::magic_name("core__walk_to_find_source_pos_info"),true,&core::core__walk_to_find_source_pos_info,"(arg &optional default-spi)");
}
NOINLINE void expose_function_687_helper() {
  expose_function(core::magic_name("core__walk_to_assign_source_pos_info"),true,&core::core__walk_to_assign_source_pos_info,"(obj top &optional stream)");
}
NOINLINE void expose_function_688_helper() {
  expose_function(core::magic_name("core__walk_to_find_source_info"),true,&core::core__walk_to_find_source_info,"(arg)");
}
NOINLINE void expose_function_689_helper() {
  expose_function(core::magic_name("core__source_pos_info_column"),true,&core::core__source_pos_info_column,"(source-pos-info)");
}
NOINLINE void expose_function_690_helper() {
  expose_function(core::magic_name("core__source_pos_info_lineno"),true,&core::core__source_pos_info_lineno,"(source-pos-info)");
}
NOINLINE void expose_function_691_helper() {
  expose_function(core::magic_name("core__source_pos_info_filepos"),true,&core::core__source_pos_info_filepos,"(source-pos-info)");
}
NOINLINE void expose_function_692_helper() {
  expose_function(core::magic_name("core__source_file_info"),true,&core::core__source_file_info,"(name &optional source-debug-namestring (source-debug-offset 0) (use-lineno t))");
}
NOINLINE void expose_function_693_helper() {
  expose_function(core::magic_name("cl__print"),true,&core::cl__print,"(obj &optional output-stream-desig)");
}
NOINLINE void expose_function_694_helper() {
  expose_function(core::magic_name("cl__prin1"),true,&core::cl__prin1,"(obj &optional output-stream-desig)");
}
NOINLINE void expose_function_695_helper() {
  expose_function(core::magic_name("cl__princ"),true,&core::cl__princ,"(obj &optional output-stream-desig)");
}
NOINLINE void expose_function_696_helper() {
  expose_function(core::magic_name("cl__pprint"),true,&core::cl__pprint,"(obj &optional stream)");
}
NOINLINE void expose_function_697_helper() {
  expose_function(core::magic_name("core__print_unreadable_object_function"),true,&core::core__print_unreadable_object_function,"(o stream type id function)");
}
NOINLINE void expose_function_698_helper() {
  expose_function(core::magic_name("cl__write"),true,&core::cl__write,"(x &key ((:stream strm) nil) (array *print-array*) (base *print-base*) ((:case cas) *print-case*) (circle *print-circle*) (escape *print-escape*) (gensym *print-gensym*) (length *print-length*) (level *print-level*) (lines *print-lines*) (miser_width *print-miser-width*) (pprint_dispatch *print-pprint-dispatch*) (pretty *print-pretty*) (radix *print-radix*) (readably *print-readably*) (right_margin *print-right-margin*))");
}
NOINLINE void expose_function_699_helper() {
  expose_function(core::magic_name("core__write_object"),true,&core::core__write_object,"(obj &optional strm)");
}
NOINLINE void expose_function_700_helper() {
  expose_function(core::magic_name("make-bignum","CORE"),true,&core::Bignum_O::make,"");
}
NOINLINE void expose_function_701_helper() {
  expose_function(core::magic_name("cl__load"),true,&core::cl__load,"(source &key (verbose *load-verbose*) (print *load-print*) (if-does-not-exist :error) (external-format :default) (search-list core::*load-search-list*))");
}
NOINLINE void expose_function_702_helper() {
  expose_function(core::magic_name("core__load_source"),true,&core::core__load_source,"(source &optional verbose print external-format)");
}
NOINLINE void expose_function_703_helper() {
  expose_function(core::magic_name("cl__imagpart"),true,&core::cl__imagpart,"(x)");
}
NOINLINE void expose_function_704_helper() {
  expose_function(core::magic_name("cl__realpart"),true,&core::cl__realpart,"(x)");
}
NOINLINE void expose_function_705_helper() {
  expose_function(core::magic_name("cl__complex"),true,&core::cl__complex,"(r &optional (i 0))");
}
NOINLINE void expose_function_706_helper() {
  expose_function(core::magic_name("cl__integer_decode_float"),true,&core::cl__integer_decode_float,"(x)");
}
NOINLINE void expose_function_707_helper() {
  expose_function(core::magic_name("cl__float_precision"),true,&core::cl__float_precision,"(Float_sp x)");
}
NOINLINE void expose_function_708_helper() {
  expose_function(core::magic_name("cl__float_digits"),true,&core::cl__float_digits,"(x)");
}
NOINLINE void expose_function_709_helper() {
  expose_function(core::magic_name("cl__float_sign"),true,&core::cl__float_sign,"(x &optional (y x yp))");
}
NOINLINE void expose_function_710_helper() {
  expose_function(core::magic_name("cl__scale_float"),true,&core::cl__scale_float,"(x y)");
}
NOINLINE void expose_function_711_helper() {
  expose_function(core::magic_name("cl__decode_float"),true,&core::cl__decode_float,"(x)");
}
NOINLINE void expose_function_712_helper() {
  expose_function(core::magic_name("cl__rem"),true,&core::cl__rem,"(x y)");
}
NOINLINE void expose_function_713_helper() {
  expose_function(core::magic_name("cl__mod"),true,&core::cl__mod,"(x y)");
}
NOINLINE void expose_function_714_helper() {
  expose_function(core::magic_name("cl__round"),true,&core::cl__round,"(x &optional y)");
}
NOINLINE void expose_function_715_helper() {
  expose_function(core::magic_name("cl__truncate"),true,&core::cl__truncate,"(x &optional y)");
}
NOINLINE void expose_function_716_helper() {
  expose_function(core::magic_name("cl__ceiling"),true,&core::cl__ceiling,"(x &optional y)");
}
NOINLINE void expose_function_717_helper() {
  expose_function(core::magic_name("cl__floor"),true,&core::cl__floor,"(x &optional y)");
}
NOINLINE void expose_function_718_helper() {
  expose_function(core::magic_name("cl__denominator"),true,&core::cl__denominator,"(x)");
}
NOINLINE void expose_function_719_helper() {
  expose_function(core::magic_name("cl__numerator"),true,&core::cl__numerator,"(x)");
}
NOINLINE void expose_function_720_helper() {
  expose_function(core::magic_name("cl__float"),true,&core::cl__float,"(x &optional y)");
}
NOINLINE void expose_function_721_helper() {
  expose_function(core::magic_name("core__integer_to_string"),true,&core::core__integer_to_string,"(buffer integer base radix decimalp)");
}
NOINLINE void expose_function_722_helper() {
  expose_function(core::magic_name("core__bignum_to_string"),true,&core::core__bignum_to_string,"(buffer x base)");
}
NOINLINE void expose_function_723_helper() {
  expose_function(core::magic_name("cl__lcm"),true,&core::cl__lcm,"(&rest args)");
}
NOINLINE void expose_function_724_helper() {
  expose_function(core::magic_name("cl__gcd"),true,&core::cl__gcd,"(&rest nums)");
}
NOINLINE void expose_function_725_helper() {
  expose_function(core::magic_name("core__num_op_atanh"),true,&core::core__num_op_atanh,"(arg)");
}
NOINLINE void expose_function_726_helper() {
  expose_function(core::magic_name("core__num_op_acosh"),true,&core::core__num_op_acosh,"(arg)");
}
NOINLINE void expose_function_727_helper() {
  expose_function(core::magic_name("core__num_op_asinh"),true,&core::core__num_op_asinh,"(arg)");
}
NOINLINE void expose_function_728_helper() {
  expose_function(core::magic_name("core__num_op_acos"),true,&core::core__num_op_acos,"(arg)");
}
NOINLINE void expose_function_729_helper() {
  expose_function(core::magic_name("core__num_op_asin"),true,&core::core__num_op_asin,"(arg)");
}
NOINLINE void expose_function_730_helper() {
  expose_function(core::magic_name("CORE:randomNumberNormal01"),true,&core::randomNumberNormal01,"");
}
NOINLINE void expose_function_731_helper() {
  expose_function(core::magic_name("CORE:randomNumber01"),true,&core::randomNumber01,"");
}
NOINLINE void expose_function_732_helper() {
  expose_function(core::magic_name("CORE:seedRandomNumberGeneratorsUsingTime"),true,&core::seedRandomNumberGeneratorsUsingTime,"");
}
NOINLINE void expose_function_733_helper() {
  expose_function(core::magic_name("CORE:seedRandomNumberGenerators"),true,&core::seedRandomNumberGenerators,"");
}
NOINLINE void expose_function_734_helper() {
  expose_function(core::magic_name("cl__get_universal_time"),true,&core::cl__get_universal_time,"");
}
NOINLINE void expose_function_735_helper() {
  expose_function(core::magic_name("CORE:general-two-arg-_EQ_"),true,&core::two_arg__EQ_,"");
}
NOINLINE void expose_function_736_helper() {
  expose_function(core::magic_name("CORE:general-two-arg-_GE_"),true,&core::two_arg__GE_,"");
}
NOINLINE void expose_function_737_helper() {
  expose_function(core::magic_name("CORE:general-two-arg-_GT_"),true,&core::two_arg__GT_,"");
}
NOINLINE void expose_function_738_helper() {
  expose_function(core::magic_name("CORE:general-two-arg-_LE_"),true,&core::two_arg__LE_,"");
}
NOINLINE void expose_function_739_helper() {
  expose_function(core::magic_name("CORE:general-two-arg-_LT_"),true,&core::two_arg__LT_,"");
}
NOINLINE void expose_function_740_helper() {
  expose_function(core::magic_name("CORE:general-two-arg-_DIVIDE_"),true,&core::contagen_div,"");
}
NOINLINE void expose_function_741_helper() {
  expose_function(core::magic_name("CORE:general-two-arg-_TIMES_"),true,&core::contagen_mul,"");
}
NOINLINE void expose_function_742_helper() {
  expose_function(core::magic_name("CORE:general-two-arg-_MINUS_"),true,&core::contagen_sub,"");
}
NOINLINE void expose_function_743_helper() {
  expose_function(core::magic_name("CORE:general-two-arg-_PLUS_"),true,&core::contagen_add,"");
}
NOINLINE void expose_function_744_helper() {
  expose_function(core::magic_name("cl__integer_length"),true,&core::cl__integer_length,"(i)");
}
NOINLINE void expose_function_745_helper() {
  expose_function(core::magic_name("core__log1p"),true,&core::core__log1p,"(arg)");
}
NOINLINE void expose_function_746_helper() {
  expose_function(core::magic_name("cl__log"),true,&core::cl__log,"(number &optional base)");
}
NOINLINE void expose_function_747_helper() {
  expose_function(core::magic_name("cl__atan"),true,&core::cl__atan,"(x &optional y)");
}
NOINLINE void expose_function_748_helper() {
  expose_function(core::magic_name("cl__expt"),true,&core::cl__expt,"(x y)");
}
NOINLINE void expose_function_749_helper() {
  expose_function(core::magic_name("cl__exp"),true,&core::cl__exp,"(x)");
}
NOINLINE void expose_function_750_helper() {
  expose_function(core::magic_name("cl__conjugate"),true,&core::cl__conjugate,"(x)");
}
NOINLINE void expose_function_751_helper() {
  expose_function(core::magic_name("cl__tanh"),true,&core::cl__tanh,"(x)");
}
NOINLINE void expose_function_752_helper() {
  expose_function(core::magic_name("cl__cosh"),true,&core::cl__cosh,"(x)");
}
NOINLINE void expose_function_753_helper() {
  expose_function(core::magic_name("cl__sinh"),true,&core::cl__sinh,"(x)");
}
NOINLINE void expose_function_754_helper() {
  expose_function(core::magic_name("cl__tan"),true,&core::cl__tan,"(x)");
}
NOINLINE void expose_function_755_helper() {
  expose_function(core::magic_name("cl__cos"),true,&core::cl__cos,"(x)");
}
NOINLINE void expose_function_756_helper() {
  expose_function(core::magic_name("cl__sin"),true,&core::cl__sin,"(x)");
}
NOINLINE void expose_function_757_helper() {
  expose_function(core::magic_name("cl__sqrt"),true,&core::cl__sqrt,"(arg)");
}
NOINLINE void expose_function_758_helper() {
  expose_function(core::magic_name("core__nan"),true,&core::core__nan,"(num)");
}
NOINLINE void expose_function_759_helper() {
  expose_function(core::magic_name("cl___EQ_"),true,&core::cl___EQ_,"(&rest args)");
}
NOINLINE void expose_function_760_helper() {
  expose_function(core::magic_name("cl___NE_"),true,&core::cl___NE_,"(&rest args)");
}
NOINLINE void expose_function_761_helper() {
  expose_function(core::magic_name("cl___GE_"),true,&core::cl___GE_,"(&rest args)");
}
NOINLINE void expose_function_762_helper() {
  expose_function(core::magic_name("cl___LE_"),true,&core::cl___LE_,"(&rest args)");
}
NOINLINE void expose_function_763_helper() {
  expose_function(core::magic_name("cl___GT_"),true,&core::cl___GT_,"(&rest args)");
}
NOINLINE void expose_function_764_helper() {
  expose_function(core::magic_name("cl___LT_"),true,&core::cl___LT_,"(&rest args)");
}
NOINLINE void expose_function_765_helper() {
  expose_function("CORE:TWO-ARG-=",true,&core::two_arg__EQ_,"");
}
NOINLINE void expose_function_766_helper() {
  expose_function("CORE:TWO-ARG->=",true,&core::two_arg__GE_,"");
}
NOINLINE void expose_function_767_helper() {
  expose_function("CORE:TWO-ARG->",true,&core::two_arg__GT_,"");
}
NOINLINE void expose_function_768_helper() {
  expose_function("CORE:TWO-ARG-<=",true,&core::two_arg__LE_,"");
}
NOINLINE void expose_function_769_helper() {
  expose_function("CORE:TWO-ARG-<",true,&core::two_arg__LT_,"");
}
NOINLINE void expose_function_770_helper() {
  expose_function(core::magic_name("cl___DIVIDE_"),true,&core::cl___DIVIDE_,"(num &rest numbers)");
}
NOINLINE void expose_function_771_helper() {
  expose_function(core::magic_name("cl___MINUS_"),true,&core::cl___MINUS_,"(num &rest numbers)");
}
NOINLINE void expose_function_772_helper() {
  expose_function(core::magic_name("cl___TIMES_"),true,&core::cl___TIMES_,"(&rest numbers)");
}
NOINLINE void expose_function_773_helper() {
  expose_function(core::magic_name("cl___PLUS_"),true,&core::cl___PLUS_,"(&rest numbers)");
}
NOINLINE void expose_function_774_helper() {
  expose_function("CORE:TWO-ARG-/",true,&core::contagen_div,"");
}
NOINLINE void expose_function_775_helper() {
  expose_function("CORE:TWO-ARG-*",true,&core::contagen_mul,"");
}
NOINLINE void expose_function_776_helper() {
  expose_function("CORE:TWO-ARG--",true,&core::contagen_sub,"");
}
NOINLINE void expose_function_777_helper() {
  expose_function("CORE:TWO-ARG-+",true,&core::contagen_add,"");
}
NOINLINE void expose_function_778_helper() {
  expose_function(core::magic_name("cl__lognor"),true,&core::cl__lognor,"(a b)");
}
NOINLINE void expose_function_779_helper() {
  expose_function(core::magic_name("cl__lognand"),true,&core::cl__lognand,"(a b)");
}
NOINLINE void expose_function_780_helper() {
  expose_function(core::magic_name("cl__lognot"),true,&core::cl__lognot,"(a)");
}
NOINLINE void expose_function_781_helper() {
  expose_function(core::magic_name("cl__logorc2"),true,&core::cl__logorc2,"(a b)");
}
NOINLINE void expose_function_782_helper() {
  expose_function(core::magic_name("cl__logorc1"),true,&core::cl__logorc1,"(a b)");
}
NOINLINE void expose_function_783_helper() {
  expose_function(core::magic_name("cl__logandc2"),true,&core::cl__logandc2,"(a b)");
}
NOINLINE void expose_function_784_helper() {
  expose_function(core::magic_name("cl__logandc1"),true,&core::cl__logandc1,"(a b)");
}
NOINLINE void expose_function_785_helper() {
  expose_function(core::magic_name("cl__logeqv"),true,&core::cl__logeqv,"(&rest integers)");
}
NOINLINE void expose_function_786_helper() {
  expose_function(core::magic_name("cl__logxor"),true,&core::cl__logxor,"(&rest integers)");
}
NOINLINE void expose_function_787_helper() {
  expose_function(core::magic_name("cl__logior"),true,&core::cl__logior,"(&rest integers)");
}
NOINLINE void expose_function_788_helper() {
  expose_function(core::magic_name("cl__logand"),true,&core::cl__logand,"(&rest integers)");
}
NOINLINE void expose_function_789_helper() {
  expose_function(core::magic_name("cl__max"),true,&core::cl__max,"(max &rest nums)");
}
NOINLINE void expose_function_790_helper() {
  expose_function(core::magic_name("cl__min"),true,&core::cl__min,"(min &rest nums)");
}
NOINLINE void expose_function_791_helper() {
  expose_function(core::magic_name("core__fixnum_number_of_bits"),true,&core::core__fixnum_number_of_bits,"");
}
NOINLINE void expose_function_792_helper() {
  expose_function(core::magic_name("core__convert_overflow_result_to_bignum"),true,&core::core__convert_overflow_result_to_bignum,"");
}
NOINLINE void expose_function_793_helper() {
  expose_function(core::magic_name("cl__zerop"),true,&core::cl__zerop,"(num)");
}
NOINLINE void expose_function_794_helper() {
  expose_function(core::magic_name("core__hash_table_force_rehash"),true,&core::core__hash_table_force_rehash,"(ht)");
}
NOINLINE void expose_function_795_helper() {
  expose_function(core::magic_name("cl__gethash"),true,&core::cl__gethash,"(key hash-table &optional default_value)");
}
NOINLINE void expose_function_796_helper() {
  expose_function(core::magic_name("cl__remhash"),true,&core::cl__remhash,"(key hashtable)");
}
NOINLINE void expose_function_797_helper() {
  expose_function(core::magic_name("core__hash_equalp"),true,&core::core__hash_equalp,"(&rest args)");
}
NOINLINE void expose_function_798_helper() {
  expose_function(core::magic_name("core__hash_equal"),true,&core::core__hash_equal,"(&rest args)");
}
NOINLINE void expose_function_799_helper() {
  expose_function(core::magic_name("core__hash_eql"),true,&core::core__hash_eql,"(&rest args)");
}
NOINLINE void expose_function_800_helper() {
  expose_function(core::magic_name("core__hash_table_entry_deleted_p"),true,&core::core__hash_table_entry_deleted_p,"(cons)");
}
NOINLINE void expose_function_801_helper() {
  expose_function(core::magic_name("cl__clrhash"),true,&core::cl__clrhash,"(hash_table)");
}
NOINLINE void expose_function_802_helper() {
  expose_function(core::magic_name("cl__maphash"),true,&core::cl__maphash,"(function_desig hash_table)");
}
NOINLINE void expose_function_803_helper() {
  expose_function(core::magic_name("core__hash_table_weakness"),true,&core::core__hash_table_weakness,"(ht)");
}
NOINLINE void expose_function_804_helper() {
  expose_function(core::magic_name("cl__make_hash_table"),true,&core::cl__make_hash_table,"(&key (test (function eql)) (size 16) (rehash-size 1.5) (rehash_threshold 1.0) weakness debug)");
}
NOINLINE void expose_function_805_helper() {
  expose_function(core::magic_name("core__debug_hash_table"),true,&core::core__debug_hash_table,"(on)");
}
NOINLINE void expose_function_806_helper() {
  expose_function(core::magic_name("core__make_small_multimap"),true,&core::core__make_small_multimap,"");
}
NOINLINE void expose_function_807_helper() {
  expose_function(core::magic_name("core__make_small_map"),true,&core::core__make_small_map,"");
}
NOINLINE void expose_function_808_helper() {
  expose_function(core::magic_name("core__print_current_ihs_frame_environment"),true,&core::core__print_current_ihs_frame_environment,"");
}
NOINLINE void expose_function_809_helper() {
  expose_function(core::magic_name("core__frame_pointers"),true,&core::core__frame_pointers,"");
}
NOINLINE void expose_function_810_helper() {
  expose_function(core::magic_name("core__clib_backtrace"),true,&core::core__clib_backtrace,"(depth)");
}
NOINLINE void expose_function_811_helper() {
  expose_function(core::magic_name("core__low_level_backtrace"),true,&core::core__low_level_backtrace,"");
}
NOINLINE void expose_function_812_helper() {
  expose_function(core::magic_name("cl__random"),true,&core::cl__random,"(olimit &optional (random-state cl:*random-state*))");
}
NOINLINE void expose_function_813_helper() {
  expose_function(core::magic_name("make-random-state","COMMON-LISP"),true,&core::RandomState_O::make,"(&optional state)");
}
NOINLINE void expose_function_814_helper() {
  expose_function(core::magic_name("core__ihs_backtrace"),true,&core::core__ihs_backtrace,"(&optional (out t) msg)");
}
NOINLINE void expose_function_815_helper() {
  expose_function(core::magic_name("core__dynamic_binding_stack_dump"),true,&core::core__dynamic_binding_stack_dump,"");
}
NOINLINE void expose_function_816_helper() {
  expose_function(core::magic_name("core__exception_stack_dump"),true,&core::core__exception_stack_dump,"");
}
NOINLINE void expose_function_817_helper() {
  expose_function(core::magic_name("core__exception_stack"),true,&core::core__exception_stack,"");
}
NOINLINE void expose_function_818_helper() {
  expose_function(core::magic_name("core__bds_val"),true,&core::core__bds_val,"(idx)");
}
NOINLINE void expose_function_819_helper() {
  expose_function(core::magic_name("core__bds_var"),true,&core::core__bds_var,"(idx)");
}
NOINLINE void expose_function_820_helper() {
  expose_function(core::magic_name("core__bds_top"),true,&core::core__bds_top,"");
}
NOINLINE void expose_function_821_helper() {
  expose_function(core::magic_name("core__set_ihs_current_frame"),true,&core::core__set_ihs_current_frame,"");
}
NOINLINE void expose_function_822_helper() {
  expose_function(core::magic_name("core__ihs_current_frame"),true,&core::core__ihs_current_frame,"");
}
NOINLINE void expose_function_823_helper() {
  expose_function(core::magic_name("core__ihs_bds"),true,&core::core__ihs_bds,"(cur)");
}
NOINLINE void expose_function_824_helper() {
  expose_function(core::magic_name("core__ihs_env"),true,&core::core__ihs_env,"(cur)");
}
NOINLINE void expose_function_825_helper() {
  expose_function(core::magic_name("core__ihs_arguments"),true,&core::core__ihs_arguments,"(arg)");
}
NOINLINE void expose_function_826_helper() {
  expose_function(core::magic_name("core__ihs_fun"),true,&core::core__ihs_fun,"(arg)");
}
NOINLINE void expose_function_827_helper() {
  expose_function(core::magic_name("core__ihs_next"),true,&core::core__ihs_next,"(cur)");
}
NOINLINE void expose_function_828_helper() {
  expose_function(core::magic_name("core__ihs_prev"),true,&core::core__ihs_prev,"(cur)");
}
NOINLINE void expose_function_829_helper() {
  expose_function(core::magic_name("core__ihs_top"),true,&core::core__ihs_top,"");
}
NOINLINE void expose_function_830_helper() {
  expose_function(core::magic_name("core__ihs_backtrace_no_args"),true,&core::core__ihs_backtrace_no_args,"");
}
NOINLINE void expose_function_831_helper() {
  expose_function(core::magic_name("core__get_invocation_history_frame_next"),true,&core::core__get_invocation_history_frame_next,"");
}
NOINLINE void expose_function_832_helper() {
  expose_function(core::magic_name("core__get_invocation_history_frame_prev"),true,&core::core__get_invocation_history_frame_prev,"");
}
NOINLINE void expose_function_833_helper() {
  expose_function(core::magic_name("core__get_invocation_history_frame"),true,&core::core__get_invocation_history_frame,"");
}
NOINLINE void expose_function_834_helper() {
  expose_function(core::magic_name("core__get_invocation_history_frame_top"),true,&core::core__get_invocation_history_frame_top,"");
}
NOINLINE void expose_function_835_helper() {
  expose_function(core::magic_name("core__get_invocation_history_frame_search"),true,&core::core__get_invocation_history_frame_search,"(idx direction)");
}
NOINLINE void expose_function_836_helper() {
  expose_function(core::magic_name("CORE:make-invocation-history-frame-iterator"),true,&core::InvocationHistoryFrameIterator_O::make,"");
}
NOINLINE void expose_function_837_helper() {
  expose_function(core::magic_name("cl__sxhash"),true,&core::cl__sxhash,"(obj)");
}
NOINLINE void expose_function_838_helper() {
  expose_function(core::magic_name("cl__type_of"),true,&core::cl__type_of,"(obj)");
}
NOINLINE void expose_function_839_helper() {
  expose_function(core::magic_name("core__type_to_symbol"),true,&core::core__type_to_symbol,"(x)");
}
NOINLINE void expose_function_840_helper() {
  expose_function(core::magic_name("cl__gensym"),true,&core::cl__gensym,"(&optional x)");
}
NOINLINE void expose_function_841_helper() {
  expose_function(core::magic_name("core__sequence_start_end"),true,&core::core__sequence_start_end,"(func sequence start end)");
}
NOINLINE void expose_function_842_helper() {
  expose_function(core::magic_name("cl__append"),true,&core::cl__append,"(&rest lists)");
}
NOINLINE void expose_function_843_helper() {
  expose_function(core::magic_name("cl__mapcan"),true,&core::cl__mapcan,"(op &rest lists)");
}
NOINLINE void expose_function_844_helper() {
  expose_function(core::magic_name("cl__mapcon"),true,&core::cl__mapcon,"(op &rest lists)");
}
NOINLINE void expose_function_845_helper() {
  expose_function(core::magic_name("core__mapappend"),true,&core::core__mapappend,"(fun &rest cargs)");
}
NOINLINE void expose_function_846_helper() {
  expose_function(core::magic_name("cl__mapl"),true,&core::cl__mapl,"(op &rest lists)");
}
NOINLINE void expose_function_847_helper() {
  expose_function(core::magic_name("cl__maplist"),true,&core::cl__maplist,"(func_desig &rest lists)");
}
NOINLINE void expose_function_848_helper() {
  expose_function(core::magic_name("cl__mapc"),true,&core::cl__mapc,"(op &rest lists)");
}
NOINLINE void expose_function_849_helper() {
  expose_function(core::magic_name("cl__mapcar"),true,&core::cl__mapcar,"(func_desig &rest lists)");
}
NOINLINE void expose_function_850_helper() {
  expose_function(core::magic_name("cl__notevery"),true,&core::cl__notevery,"(predicate &rest sequences)");
}
NOINLINE void expose_function_851_helper() {
  expose_function(core::magic_name("cl__notany"),true,&core::cl__notany,"(predicate &rest sequences)");
}
NOINLINE void expose_function_852_helper() {
  expose_function(core::magic_name("cl__some"),true,&core::cl__some,"(predicate &rest sequences)");
}
NOINLINE void expose_function_853_helper() {
  expose_function(core::magic_name("cl__every"),true,&core::cl__every,"(predicate &rest sequences)");
}
NOINLINE void expose_function_854_helper() {
  expose_function(core::magic_name("cl__read_preserving_whitespace"),true,&core::cl__read_preserving_whitespace,"(&optional input-stream-designator (eof-error-p t) eof-value recursive-p)");
}
NOINLINE void expose_function_855_helper() {
  expose_function(core::magic_name("cl__read"),true,&core::cl__read,"(&optional input-stream-designator (eof-error-p t) eof-value recursive-p)");
}
NOINLINE void expose_function_856_helper() {
  expose_function(core::magic_name("cl__read_delimited_list"),true,&core::cl__read_delimited_list,"(char &optional input-stream-designator recursive-p)");
}
NOINLINE void expose_function_857_helper() {
  expose_function(core::magic_name("cl__fmakunbound"),true,&core::cl__fmakunbound,"(function-name)");
}
NOINLINE void expose_function_858_helper() {
  expose_function(core::magic_name("cl__fboundp"),true,&core::cl__fboundp,"(function-name)");
}
NOINLINE void expose_function_859_helper() {
  expose_function(core::magic_name("cl__fdefinition"),true,&core::cl__fdefinition,"(function-name)");
}
NOINLINE void expose_function_860_helper() {
  expose_function(core::magic_name("core__STARfset"),true,&core::core__STARfset,"(function-name fn &optional macro)");
}
NOINLINE void expose_function_861_helper() {
  expose_function(core::magic_name("cl__class_of"),true,&core::cl__class_of,"(obj)");
}
NOINLINE void expose_function_862_helper() {
  expose_function(core::magic_name("cl__null"),true,&core::cl__null,"(obj)");
}
NOINLINE void expose_function_863_helper() {
  expose_function(core::magic_name("core__macroexpand_default"),true,&core::core__macroexpand_default,"(macro_function form macro_env)");
}
NOINLINE void expose_function_864_helper() {
  expose_function(core::magic_name("cl__identity"),true,&core::cl__identity,"(arg)");
}
NOINLINE void expose_function_865_helper() {
  expose_function(core::magic_name("cl__constantp"),true,&core::cl__constantp,"(obj &optional env)");
}
NOINLINE void expose_function_866_helper() {
  expose_function(core::magic_name("core__gdb_inspect"),true,&core::core__gdb_inspect,"(msg o)");
}
NOINLINE void expose_function_867_helper() {
  expose_function(core::magic_name("core__trap_execution"),true,&core::core__trap_execution,"(&optional msg)");
}
NOINLINE void expose_function_868_helper() {
  expose_function(core::magic_name("core__cxx_lambda_list_handler_create_bindings_calls"),true,&core::core__cxx_lambda_list_handler_create_bindings_calls,"");
}
NOINLINE void expose_function_869_helper() {
  expose_function(core::magic_name("core__gdb"),true,&core::core__gdb,"(&optional msg)");
}
NOINLINE void expose_function_870_helper() {
  expose_function(core::magic_name("core__break"),true,&core::core__break,"(&optional fmt-control &rest args)");
}
NOINLINE void expose_function_871_helper() {
  expose_function(core::magic_name("cl__ash"),true,&core::cl__ash,"(integer count)");
}
NOINLINE void expose_function_872_helper() {
  expose_function(core::magic_name("cl__special_operator_p"),true,&core::cl__special_operator_p,"(symbol)");
}
NOINLINE void expose_function_873_helper() {
  expose_function(core::magic_name("cl__macro_function"),true,&core::cl__macro_function,"(symbol &optional env)");
}
NOINLINE void expose_function_874_helper() {
  expose_function(core::magic_name("core__separate_pair_list"),true,&core::core__separate_pair_list,"(listOfPairs)");
}
NOINLINE void expose_function_875_helper() {
  expose_function(core::magic_name("core__test_memory_error"),true,&core::core__test_memory_error,"");
}
NOINLINE void expose_function_876_helper() {
  expose_function(core::magic_name("core__valid_function_name_p"),true,&core::core__valid_function_name_p,"(arg)");
}
NOINLINE void expose_function_877_helper() {
  expose_function(core::magic_name("core__function_block_name"),true,&core::core__function_block_name,"(functionName)");
}
NOINLINE void expose_function_878_helper() {
  expose_function(core::magic_name("cl__values_list"),true,&core::cl__values_list,"(list)");
}
NOINLINE void expose_function_879_helper() {
  expose_function(core::magic_name("core__values_testing"),true,&core::core__values_testing,"(&rest args)");
}
NOINLINE void expose_function_880_helper() {
  expose_function(core::magic_name("cl__values"),true,&core::cl__values,"(&rest args)");
}
NOINLINE void expose_function_881_helper() {
  expose_function(core::magic_name("core__smart_pointer_details"),true,&core::core__smart_pointer_details,"");
}
NOINLINE void expose_function_882_helper() {
  expose_function(core::magic_name("core__unbound"),true,&core::core__unbound,"");
}
NOINLINE void expose_function_883_helper() {
  expose_function(core::magic_name("core__substitute"),true,&core::core__substitute,"(arg)");
}
NOINLINE void expose_function_884_helper() {
  expose_function(core::magic_name("core__is_true"),true,&core::core__is_true,"(arg)");
}
NOINLINE void expose_function_885_helper() {
  expose_function(core::magic_name("core__pointer"),true,&core::core__pointer,"(arg)");
}
NOINLINE void expose_function_886_helper() {
  expose_function(core::magic_name("core__setenv"),true,&core::core__setenv,"(name arg overwrite)");
}
NOINLINE void expose_function_887_helper() {
  expose_function(core::magic_name("core__describe_cxx_object"),true,&core::core__describe_cxx_object,"(name &optional stream)");
}
NOINLINE void expose_function_888_helper() {
  expose_function(core::magic_name("ext__llvm_version"),true,&core::ext__llvm_version,"");
}
NOINLINE void expose_function_889_helper() {
  expose_function(core::magic_name("ext__getenv"),true,&core::ext__getenv,"(arg)");
}
NOINLINE void expose_function_890_helper() {
  expose_function(core::magic_name("ext__vfork_execvp"),true,&core::ext__vfork_execvp,"(call-and-arguments)");
}
NOINLINE void expose_function_891_helper() {
  expose_function(core::magic_name("ext__system"),true,&core::ext__system,"(cmd)");
}
NOINLINE void expose_function_892_helper() {
  expose_function(core::magic_name("ext__setenv"),true,&core::ext__setenv,"(name value)");
}
NOINLINE void expose_function_893_helper() {
  expose_function(core::magic_name("core__dump_tagged_fixnum"),true,&core::core__dump_tagged_fixnum,"(arg)");
}
NOINLINE void expose_function_894_helper() {
  expose_function(core::magic_name("core__from_tagged_fixnum"),true,&core::core__from_tagged_fixnum,"(val)");
}
NOINLINE void expose_function_895_helper() {
  expose_function(core::magic_name("core__to_tagged_fixnum"),true,&core::core__to_tagged_fixnum,"(arg)");
}
NOINLINE void expose_function_896_helper() {
  expose_function(core::magic_name("core__all_registered_class_names"),true,&core::core__all_registered_class_names,"");
}
NOINLINE void expose_function_897_helper() {
  expose_function(core::magic_name("core__incomplete_next_higher_power_of_2"),true,&core::core__incomplete_next_higher_power_of_2,"(arg)");
}
NOINLINE void expose_function_898_helper() {
  expose_function(core::magic_name("core__dump_address_of"),true,&core::core__dump_address_of,"(arg)");
}
NOINLINE void expose_function_899_helper() {
  expose_function(core::magic_name("cl__set"),true,&core::cl__set,"(sym value)");
}
NOINLINE void expose_function_900_helper() {
  expose_function(core::magic_name("core__argv"),true,&core::core__argv,"(idx)");
}
NOINLINE void expose_function_901_helper() {
  expose_function(core::magic_name("core__argc"),true,&core::core__argc,"");
}
NOINLINE void expose_function_902_helper() {
  expose_function(core::magic_name("cl__machine_instance"),true,&core::cl__machine_instance,"");
}
NOINLINE void expose_function_903_helper() {
  expose_function(core::magic_name("cl__machine_version"),true,&core::cl__machine_version,"");
}
NOINLINE void expose_function_904_helper() {
  expose_function(core::magic_name("cl__machine_type"),true,&core::cl__machine_type,"");
}
NOINLINE void expose_function_905_helper() {
  expose_function(core::magic_name("cl__software_version"),true,&core::cl__software_version,"");
}
NOINLINE void expose_function_906_helper() {
  expose_function(core::magic_name("cl__software_type"),true,&core::cl__software_type,"");
}
NOINLINE void expose_function_907_helper() {
  expose_function(core::magic_name("core__create_tagged_immediate_value_or_nil"),true,&core::core__create_tagged_immediate_value_or_nil,"(obj)");
}
NOINLINE void expose_function_908_helper() {
  expose_function(core::magic_name("core__lisp_implementation_id"),true,&core::core__lisp_implementation_id,"");
}
NOINLINE void expose_function_909_helper() {
  expose_function(core::magic_name("core__single_dispatch_method_cache_status"),true,&core::core__single_dispatch_method_cache_status,"");
}
NOINLINE void expose_function_910_helper() {
  expose_function(core::magic_name("core__slot_cache_status"),true,&core::core__slot_cache_status,"");
}
NOINLINE void expose_function_911_helper() {
  expose_function(core::magic_name("core__method_cache_status"),true,&core::core__method_cache_status,"");
}
NOINLINE void expose_function_912_helper() {
  expose_function(core::magic_name("core__single_dispatch_method_cache_resize"),true,&core::core__single_dispatch_method_cache_resize,"(pow)");
}
NOINLINE void expose_function_913_helper() {
  expose_function(core::magic_name("core__slot_cache_resize"),true,&core::core__slot_cache_resize,"(pow)");
}
NOINLINE void expose_function_914_helper() {
  expose_function(core::magic_name("core__method_cache_resize"),true,&core::core__method_cache_resize,"(pow)");
}
NOINLINE void expose_function_915_helper() {
  expose_function(core::magic_name("cl__lisp_implementation_version"),true,&core::cl__lisp_implementation_version,"");
}
NOINLINE void expose_function_916_helper() {
  expose_function(core::magic_name("cl__lisp_implementation_type"),true,&core::cl__lisp_implementation_type,"");
}
NOINLINE void expose_function_917_helper() {
  expose_function(core::magic_name("cl__sleep"),true,&core::cl__sleep,"(seconds)");
}
NOINLINE void expose_function_918_helper() {
  expose_function(core::magic_name("core__clear_generic_function_dispatch_cache"),true,&core::core__clear_generic_function_dispatch_cache,"");
}
NOINLINE void expose_function_919_helper() {
  expose_function(core::magic_name("CORE:allocateRawInstance"),true,&core::Instance_O::allocateRawInstance,"");
}
NOINLINE void expose_function_920_helper() {
  expose_function(core::magic_name("core__copy_instance"),true,&core::core__copy_instance,"(obj)");
}
NOINLINE void expose_function_921_helper() {
  expose_function(core::magic_name("core__instance_class_set"),true,&core::core__instance_class_set,"(instance func)");
}
NOINLINE void expose_function_922_helper() {
  expose_function(core::magic_name("clos__setFuncallableInstanceFunction"),true,&core::clos__setFuncallableInstanceFunction,"");
}
NOINLINE void expose_function_923_helper() {
  expose_function(core::magic_name("core__clear_gfun_hash"),true,&core::core__clear_gfun_hash,"(what)");
}
NOINLINE void expose_function_924_helper() {
  expose_function(core::magic_name("core__maybe_expand_generic_function_arguments"),true,&core::core__maybe_expand_generic_function_arguments,"(args)");
}
NOINLINE void expose_function_925_helper() {
  expose_function(core::magic_name("core__progv_function"),true,&core::core__progv_function,"(symbols values func)");
}
NOINLINE void expose_function_926_helper() {
  expose_function(core::magic_name("core__throw_function"),true,&core::core__throw_function,"(tag result)");
}
NOINLINE void expose_function_927_helper() {
  expose_function(core::magic_name("core__catch_function"),true,&core::core__catch_function,"(tag func)");
}
NOINLINE void expose_function_928_helper() {
  expose_function(core::magic_name("core__multiple_value_prog1_function"),true,&core::core__multiple_value_prog1_function,"(func1 func2)");
}
NOINLINE void expose_function_929_helper() {
  expose_function(core::magic_name("core__multiple_value_funcall"),true,&core::core__multiple_value_funcall,"(function-designator &rest functions)");
}
NOINLINE void expose_function_930_helper() {
  expose_function(core::magic_name("core__funwind_protect"),true,&core::core__funwind_protect,"(protected-fn cleanup-fn)");
}
NOINLINE void expose_function_931_helper() {
  expose_function(core::magic_name("core__call_with_variable_bound"),true,&core::core__call_with_variable_bound,"(sym val thunk)");
}
NOINLINE void expose_function_932_helper() {
  expose_function(core::magic_name("compiler__implicit_compile_hook_default"),true,&core::compiler__implicit_compile_hook_default,"(form &optional environment)");
}
NOINLINE void expose_function_933_helper() {
  expose_function(core::magic_name("core__dladdr"),true,&core::core__dladdr,"(addr)");
}
NOINLINE void expose_function_934_helper() {
  expose_function(core::magic_name("core__call_dl_main_function"),true,&core::core__call_dl_main_function,"(addr)");
}
NOINLINE void expose_function_935_helper() {
  expose_function(core::magic_name("core__dlsym"),true,&core::core__dlsym,"(handle name)");
}
NOINLINE void expose_function_936_helper() {
  expose_function(core::magic_name("core__dlopen"),true,&core::core__dlopen,"(pathDesig)");
}
NOINLINE void expose_function_937_helper() {
  expose_function(core::magic_name("core__load_bundle"),true,&core::core__load_bundle,"(name &optional verbose print external-format)");
}
NOINLINE void expose_function_938_helper() {
  expose_function(core::magic_name("core__startup_image_pathname"),true,&core::core__startup_image_pathname,"");
}
NOINLINE void expose_function_939_helper() {
  expose_function(core::magic_name("core__mangle_name"),true,&core::core__mangle_name,"(object &optional is_function)");
}
NOINLINE void expose_function_940_helper() {
  expose_function(core::magic_name("core__cxx_fibn"),true,&core::core__cxx_fibn,"(reps num)");
}
NOINLINE void expose_function_941_helper() {
  expose_function(core::magic_name("core__test_tagged_cast"),true,&core::core__test_tagged_cast,"(pow2)");
}
NOINLINE void expose_function_942_helper() {
  expose_function(core::magic_name("core__help_booting"),true,&core::core__help_booting,"");
}
NOINLINE void expose_function_943_helper() {
  expose_function(core::magic_name("core__eval_with_env_default"),true,&core::core__eval_with_env_default,"(form &optional env)");
}
NOINLINE void expose_function_944_helper() {
  expose_function(core::magic_name("core__evaluate_verbosity"),true,&core::core__evaluate_verbosity,"(arg)");
}
NOINLINE void expose_function_945_helper() {
  expose_function(core::magic_name("core__evaluate_depth"),true,&core::core__evaluate_depth,"");
}
NOINLINE void expose_function_946_helper() {
  expose_function(core::magic_name("core__classify_let_variables_and_declares"),true,&core::core__classify_let_variables_and_declares,"(variables declared-specials)");
}
NOINLINE void expose_function_947_helper() {
  expose_function(core::magic_name("core__lookup_symbol_macro"),true,&core::core__lookup_symbol_macro,"(symbol &optional env)");
}
NOINLINE void expose_function_948_helper() {
  expose_function(core::magic_name("core__extract_lambda_name"),true,&core::core__extract_lambda_name,"(form &optional default)");
}
NOINLINE void expose_function_949_helper() {
  expose_function(core::magic_name("core__extract_lambda_name_from_declares"),true,&core::core__extract_lambda_name_from_declares,"(declare-list &optional default)");
}
NOINLINE void expose_function_950_helper() {
  expose_function(core::magic_name("core__process_declarations"),true,&core::core__process_declarations,"(body expectDocString)");
}
NOINLINE void expose_function_951_helper() {
  expose_function(core::magic_name("core__coerce_to_function"),true,&core::core__coerce_to_function,"(arg)");
}
NOINLINE void expose_function_952_helper() {
  expose_function(core::magic_name("cl__funcall"),true,&core::cl__funcall,"(function_desig &va-rest args)");
}
NOINLINE void expose_function_953_helper() {
  expose_function(core::magic_name("cl__eval"),true,&core::cl__eval,"(form)");
}
NOINLINE void expose_function_954_helper() {
  expose_function(core::magic_name("cl__apply"),true,&core::cl__apply,"(head &va-rest args)");
}
NOINLINE void expose_function_955_helper() {
  expose_function(core::magic_name("core__compile_form_and_eval_with_env"),true,&core::core__compile_form_and_eval_with_env,"(form &optional env stepping compiler-env-p (execute t))");
}
NOINLINE void expose_function_956_helper() {
  expose_function(core::magic_name("CORE:make-weak-pointer"),true,&core::WeakPointer_O::make,"");
}
NOINLINE void expose_function_957_helper() {
  expose_function(core::magic_name("core__weak_rehash"),true,&core::core__weak_rehash,"(ht &optional sz)");
}
NOINLINE void expose_function_958_helper() {
  expose_function(core::magic_name("core__weak_splat"),true,&core::core__weak_splat,"(ht idx)");
}
NOINLINE void expose_function_959_helper() {
  expose_function(core::magic_name("core__weak_clrhash"),true,&core::core__weak_clrhash,"(ht)");
}
NOINLINE void expose_function_960_helper() {
  expose_function(core::magic_name("core__weak_remhash"),true,&core::core__weak_remhash,"(ht key)");
}
NOINLINE void expose_function_961_helper() {
  expose_function(core::magic_name("core__weak_setf_gethash"),true,&core::core__weak_setf_gethash,"(ht key value)");
}
NOINLINE void expose_function_962_helper() {
  expose_function(core::magic_name("core__weak_gethash"),true,&core::core__weak_gethash,"(key hash-table &optional default-value)");
}
NOINLINE void expose_function_963_helper() {
  expose_function(core::magic_name("core__make_weak_key_hash_table"),true,&core::core__make_weak_key_hash_table,"(&optional (size 16))");
}
NOINLINE void expose_function_964_helper() {
  expose_function(core::magic_name("CORE:make-weak-key-mapping"),true,&core::WeakKeyMapping_O::make,"");
}
NOINLINE void expose_function_965_helper() {
  expose_function(core::magic_name("core__make_vector_displaced"),true,&core::core__make_vector_displaced,"(dim element-type displaced-to displaced-offset)");
}
NOINLINE void expose_function_966_helper() {
  expose_function(core::magic_name("core__test_add"),true,&core::core__test_add,"(x y)");
}
NOINLINE void expose_function_967_helper() {
  expose_function(core::magic_name("cl__listp"),true,&core::cl__listp,"");
}
NOINLINE void expose_function_968_helper() {
  expose_function(core::magic_name("cl__characterp"),true,&core::cl__characterp,"");
}
NOINLINE void expose_function_969_helper() {
  expose_function(core::magic_name("cl__consp"),true,&core::cl__consp,"");
}
NOINLINE void expose_function_970_helper() {
  expose_function(core::magic_name("cl__atom"),true,&core::cl__atom,"(arg)");
}
NOINLINE void expose_function_971_helper() {
  expose_function(core::magic_name("cl__endp"),true,&core::cl__endp,"(arg)");
}
NOINLINE void expose_function_972_helper() {
  expose_function(core::magic_name("core__fixnump"),true,&core::core__fixnump,"");
}
NOINLINE void expose_function_973_helper() {
  expose_function(core::magic_name("cl__symbolp"),true,&core::cl__symbolp,"");
}
NOINLINE void expose_function_974_helper() {
  expose_function(core::magic_name("CORE:negate"),true,&core::clasp_negate,"");
}
NOINLINE void expose_function_975_helper() {
  expose_function(core::magic_name("CORE:oneMinus"),true,&core::clasp_one_minus,"");
}
NOINLINE void expose_function_976_helper() {
  expose_function(core::magic_name("CORE:onePlus"),true,&core::clasp_one_plus,"");
}
NOINLINE void expose_function_977_helper() {
  expose_function(core::magic_name("signum","COMMON-LISP"),true,&core::clasp_signum,"");
}
NOINLINE void expose_function_978_helper() {
  expose_function(core::magic_name("abs","COMMON-LISP"),true,&core::clasp_abs,"");
}
NOINLINE void expose_function_979_helper() {
  expose_function(core::magic_name("oddp","COMMON-LISP"),true,&core::clasp_oddp,"");
}
NOINLINE void expose_function_980_helper() {
  expose_function(core::magic_name("evenp","COMMON-LISP"),true,&core::clasp_evenp,"");
}
NOINLINE void expose_function_981_helper() {
  expose_function(core::magic_name("minusp","COMMON-LISP"),true,&core::clasp_minusp,"");
}
NOINLINE void expose_function_982_helper() {
  expose_function(core::magic_name("plusp","COMMON-LISP"),true,&core::clasp_plusp,"");
}
NOINLINE void expose_function_983_helper() {
  expose_function(core::magic_name("Tenth","COMMON-LISP"),true,&core::oTenth,"");
}
NOINLINE void expose_function_984_helper() {
  expose_function(core::magic_name("Ninth","COMMON-LISP"),true,&core::oNinth,"");
}
NOINLINE void expose_function_985_helper() {
  expose_function(core::magic_name("Eighth","COMMON-LISP"),true,&core::oEighth,"");
}
NOINLINE void expose_function_986_helper() {
  expose_function(core::magic_name("Seventh","COMMON-LISP"),true,&core::oSeventh,"");
}
NOINLINE void expose_function_987_helper() {
  expose_function(core::magic_name("Sixth","COMMON-LISP"),true,&core::oSixth,"");
}
NOINLINE void expose_function_988_helper() {
  expose_function(core::magic_name("Fifth","COMMON-LISP"),true,&core::oFifth,"");
}
NOINLINE void expose_function_989_helper() {
  expose_function(core::magic_name("Fourth","COMMON-LISP"),true,&core::oFourth,"");
}
NOINLINE void expose_function_990_helper() {
  expose_function(core::magic_name("Third","COMMON-LISP"),true,&core::oThird,"");
}
NOINLINE void expose_function_991_helper() {
  expose_function(core::magic_name("Second","COMMON-LISP"),true,&core::oSecond,"");
}
NOINLINE void expose_function_992_helper() {
  expose_function(core::magic_name("First","COMMON-LISP"),true,&core::oFirst,"");
}
NOINLINE void expose_function_993_helper() {
  expose_function(core::magic_name("cddddr","COMMON-LISP"),true,&core::oCddddr,"");
}
NOINLINE void expose_function_994_helper() {
  expose_function(core::magic_name("cddadr","COMMON-LISP"),true,&core::oCddadr,"");
}
NOINLINE void expose_function_995_helper() {
  expose_function(core::magic_name("cdaddr","COMMON-LISP"),true,&core::oCdaddr,"");
}
NOINLINE void expose_function_996_helper() {
  expose_function(core::magic_name("cdaadr","COMMON-LISP"),true,&core::oCdaadr,"");
}
NOINLINE void expose_function_997_helper() {
  expose_function(core::magic_name("cadddr","COMMON-LISP"),true,&core::oCadddr,"");
}
NOINLINE void expose_function_998_helper() {
  expose_function(core::magic_name("cadadr","COMMON-LISP"),true,&core::oCadadr,"");
}
NOINLINE void expose_function_999_helper() {
  expose_function(core::magic_name("caaddr","COMMON-LISP"),true,&core::oCaaddr,"");
}
NOINLINE void expose_function_1000_helper() {
  expose_function(core::magic_name("caaadr","COMMON-LISP"),true,&core::oCaaadr,"");
}
NOINLINE void expose_function_1001_helper() {
  expose_function(core::magic_name("cdddar","COMMON-LISP"),true,&core::oCdddar,"");
}
NOINLINE void expose_function_1002_helper() {
  expose_function(core::magic_name("cddaar","COMMON-LISP"),true,&core::oCddaar,"");
}
NOINLINE void expose_function_1003_helper() {
  expose_function(core::magic_name("cdadar","COMMON-LISP"),true,&core::oCdadar,"");
}
NOINLINE void expose_function_1004_helper() {
  expose_function(core::magic_name("cdaaar","COMMON-LISP"),true,&core::oCdaaar,"");
}
NOINLINE void expose_function_1005_helper() {
  expose_function(core::magic_name("caddar","COMMON-LISP"),true,&core::oCaddar,"");
}
NOINLINE void expose_function_1006_helper() {
  expose_function(core::magic_name("cadaar","COMMON-LISP"),true,&core::oCadaar,"");
}
NOINLINE void expose_function_1007_helper() {
  expose_function(core::magic_name("caadar","COMMON-LISP"),true,&core::oCaadar,"");
}
NOINLINE void expose_function_1008_helper() {
  expose_function(core::magic_name("caaaar","COMMON-LISP"),true,&core::oCaaaar,"");
}
NOINLINE void expose_function_1009_helper() {
  expose_function(core::magic_name("cdddr","COMMON-LISP"),true,&core::oCdddr,"");
}
NOINLINE void expose_function_1010_helper() {
  expose_function(core::magic_name("cddar","COMMON-LISP"),true,&core::oCddar,"");
}
NOINLINE void expose_function_1011_helper() {
  expose_function(core::magic_name("cdadr","COMMON-LISP"),true,&core::oCdadr,"");
}
NOINLINE void expose_function_1012_helper() {
  expose_function(core::magic_name("cdaar","COMMON-LISP"),true,&core::oCdaar,"");
}
NOINLINE void expose_function_1013_helper() {
  expose_function(core::magic_name("caddr","COMMON-LISP"),true,&core::oCaddr,"");
}
NOINLINE void expose_function_1014_helper() {
  expose_function(core::magic_name("cadar","COMMON-LISP"),true,&core::oCadar,"");
}
NOINLINE void expose_function_1015_helper() {
  expose_function(core::magic_name("caadr","COMMON-LISP"),true,&core::oCaadr,"");
}
NOINLINE void expose_function_1016_helper() {
  expose_function(core::magic_name("caaar","COMMON-LISP"),true,&core::oCaaar,"");
}
NOINLINE void expose_function_1017_helper() {
  expose_function(core::magic_name("cddr","COMMON-LISP"),true,&core::oCddr,"");
}
NOINLINE void expose_function_1018_helper() {
  expose_function(core::magic_name("cdar","COMMON-LISP"),true,&core::oCdar,"");
}
NOINLINE void expose_function_1019_helper() {
  expose_function(core::magic_name("cadr","COMMON-LISP"),true,&core::oCadr,"");
}
NOINLINE void expose_function_1020_helper() {
  expose_function(core::magic_name("caar","COMMON-LISP"),true,&core::oCaar,"");
}
NOINLINE void expose_function_1021_helper() {
  expose_function(core::magic_name("cl__rest"),true,&core::cl__rest,"");
}
NOINLINE void expose_function_1022_helper() {
  expose_function(core::magic_name("cdr","COMMON-LISP"),true,&core::oCdr,"");
}
NOINLINE void expose_function_1023_helper() {
  expose_function(core::magic_name("car","COMMON-LISP"),true,&core::oCar,"");
}
NOINLINE void expose_function_1024_helper() {
  expose_function(core::magic_name("cl__equalp"),true,&core::cl__equalp,"(x y)");
}
NOINLINE void expose_function_1025_helper() {
  expose_function(core::magic_name("cl__equal"),true,&core::cl__equal,"(x y)");
}
NOINLINE void expose_function_1026_helper() {
  expose_function(core::magic_name("cl__eql"),true,&core::cl__eql,"(x y)");
}
NOINLINE void expose_function_1027_helper() {
  expose_function(core::magic_name("cl__eq"),true,&core::cl__eq,"(x y)");
}
NOINLINE void expose_function_1028_helper() {
  expose_function(core::magic_name("ext__maybeQuote"),true,&ext::ext__maybeQuote,"");
}
NOINLINE void expose_function_1029_helper() {
  expose_function(core::magic_name("UNITS:make-quantity"),true,&units::Quantity_O::make,"");
}
NOINLINE void expose_function_1030_helper() {
  expose_function(core::magic_name("UNITS:make-unit"),true,&units::Unit_O::make,"");
}
NOINLINE void expose_function_1031_helper() {
  expose_function(core::magic_name("chem__makeIndexedObjectBag"),true,&adapt::chem__makeIndexedObjectBag,"");
}
NOINLINE void expose_function_1032_helper() {
  expose_function(core::magic_name("GEOM:make-coordinate-array"),true,&geom::CoordinateArray_O::make,"(fnsize &optional vals)");
}
NOINLINE void expose_function_1033_helper() {
  expose_function(core::magic_name("GEOM:make-omatrix"),true,&geom::OMatrix_O::make,"");
}
NOINLINE void expose_function_1034_helper() {
  expose_function(core::magic_name("GEOM:make_ovector3"),true,&geom::OVector3_O::make,"");
}
NOINLINE void expose_function_1035_helper() {
  expose_function(core::magic_name("chem__oligomer"),true,&chem::chem__oligomer,"");
}
NOINLINE void expose_function_1036_helper() {
  expose_function(core::magic_name("chem__calculate_point"),true,&chem::chem__calculate_point,"");
}
NOINLINE void expose_function_1037_helper() {
  expose_function(core::magic_name("chem__link"),true,&chem::chem__link,"");
}
NOINLINE void expose_function_1038_helper() {
  expose_function(core::magic_name("chem__monomer"),true,&chem::chem__monomer,"(monomerId groupName &optional monomerAliases comment)");
}
NOINLINE void expose_function_1039_helper() {
  expose_function(core::magic_name("chem__atom_pos"),true,&chem::chem__atom_pos,"(&rest args)");
}
NOINLINE void expose_function_1040_helper() {
  expose_function(core::magic_name("chem__find_residue"),true,&chem::chem__find_residue,"(&rest args)");
}
NOINLINE void expose_function_1041_helper() {
  expose_function(core::magic_name("chem__save_mol2"),true,&chem::chem__save_mol2,"");
}
NOINLINE void expose_function_1042_helper() {
  expose_function(core::magic_name("chem__load_mol2"),true,&chem::chem__load_mol2,"");
}
NOINLINE void expose_function_1043_helper() {
  expose_function(core::magic_name("chem__make_coordinate_array_from_atom_list"),true,&chem::chem__make_coordinate_array_from_atom_list,"");
}
NOINLINE void expose_function_1044_helper() {
  expose_function(core::magic_name("chem__save_archive_with_auto_set_cando_database"),true,&chem::chem__save_archive_with_auto_set_cando_database,"");
}
NOINLINE void expose_function_1045_helper() {
  expose_function(core::magic_name("chem__create_virtual_atom"),true,&chem::chem__create_virtual_atom,"");
}
NOINLINE void expose_function_1046_helper() {
  expose_function(core::magic_name("chem__bundle_database_path"),true,&chem::chem__bundle_database_path,"");
}
NOINLINE void expose_function_1047_helper() {
  expose_function(core::magic_name("chem__standard_database"),true,&chem::chem__standard_database,"(pathDesig &optional loadSeed (overbosity 0))");
}
NOINLINE void expose_function_1048_helper() {
  expose_function(core::magic_name("chem__database"),true,&chem::chem__database,"(fileName &optional (verbosity 0))");
}
NOINLINE void expose_function_1049_helper() {
  expose_function(core::magic_name("CHEM:make-chem-draw"),true,&chem::ChemDraw_O::make,"");
}
NOINLINE void expose_function_1050_helper() {
  expose_function(core::magic_name("chem__make_read_amber_parameters"),true,&chem::chem__make_read_amber_parameters,"");
}
NOINLINE void expose_function_1051_helper() {
  expose_function(core::magic_name("chem__set_monomer_pack"),true,&chem::chem__set_monomer_pack,"(pack-name parts &optional atom-aliases)");
}
NOINLINE void expose_function_1052_helper() {
  expose_function(core::magic_name("chem__extend_aliases"),true,&chem::chem__extend_aliases,"");
}
NOINLINE void expose_function_1053_helper() {
  expose_function(core::magic_name("chem__define_monomer_pack"),true,&chem::chem__define_monomer_pack,"(pack-name parts &optional atom-aliases)");
}
NOINLINE void expose_function_1054_helper() {
  expose_function(core::magic_name("chem__parseChimeraAtomSpecs"),true,&chem::chem__parseChimeraAtomSpecs,"");
}
NOINLINE void expose_function_1055_helper() {
  expose_function(core::magic_name("chem__alias"),true,&chem::chem__alias,"");
}
NOINLINE void expose_function_1056_helper() {
  expose_function(core::magic_name("CHEM:bonds"),true,&chem::IterateBonds_O::make,"");
}
NOINLINE void expose_function_1057_helper() {
  expose_function(core::magic_name("CHEM:residues"),true,&chem::IterateResidues_O::createForMatter,"");
}
NOINLINE void expose_function_1058_helper() {
  expose_function(core::magic_name("chem__create_for_matter"),true,&chem::chem__create_for_matter,"");
}
NOINLINE void expose_function_1059_helper() {
  expose_function(core::magic_name("CHEM:ring-bonds"),true,&chem::RingFinder_O::ringBonds,"");
}
NOINLINE void expose_function_1060_helper() {
  expose_function(core::magic_name("CHEM:identify-rings"),true,&chem::RingFinder_O::identifyRings,"");
}
NOINLINE void expose_function_1061_helper() {
  expose_function(core::magic_name("CHEM:save-pdb"),true,&chem::PdbWriter_O::savePdb,"");
}
NOINLINE void expose_function_1062_helper() {
  expose_function(core::magic_name("chem__load_pdb"),true,&chem::chem__load_pdb,"");
}
NOINLINE void expose_function_1063_helper() {
  expose_function("CHEM:MAKE-SPANNING-LOOP",true,&chem::SpanningLoop_O::make,"");
}
NOINLINE void expose_function_1064_helper() {
  expose_function(core::magic_name("chem__set_oligomer"),true,&chem::chem__set_oligomer,"");
}
NOINLINE void expose_function_1065_helper() {
  expose_function(core::magic_name("chem__oligomer_sequence"),true,&chem::chem__oligomer_sequence,"");
}
NOINLINE void expose_function_1066_helper() {
  expose_function(core::magic_name("chem__atomsInMatterAsObjectSet"),true,&chem::chem__atomsInMatterAsObjectSet,"");
}
NOINLINE void expose_function_1067_helper() {
  expose_function(core::magic_name("chem__atomsWithinSphereAsObjectSet"),true,&chem::chem__atomsWithinSphereAsObjectSet,"");
}
NOINLINE void expose_function_1068_helper() {
  expose_function(core::magic_name("chem__connectAtomsInMatterInCovalentContact"),true,&chem::chem__connectAtomsInMatterInCovalentContact,"");
}
NOINLINE void expose_function_1069_helper() {
  expose_function(core::magic_name("CHEM:otherPlugName"),true,&chem::DirectionalCoupling_O::otherPlugName,"");
}
NOINLINE void expose_function_1070_helper() {
  expose_function(core::magic_name("chem__is_in_plug_name"),true,&chem::chem__is_in_plug_name,"");
}
NOINLINE void expose_function_1071_helper() {
  expose_function(core::magic_name("CHEM:outPlugName"),true,&chem::DirectionalCoupling_O::outPlugName,"");
}
NOINLINE void expose_function_1072_helper() {
  expose_function(core::magic_name("CHEM:inPlugName"),true,&chem::DirectionalCoupling_O::inPlugName,"");
}
NOINLINE void expose_function_1073_helper() {
  expose_function(core::magic_name("CHEM:couplingName"),true,&chem::DirectionalCoupling_O::couplingName,"");
}
NOINLINE void expose_function_1074_helper() {
  expose_function(core::magic_name("CHEM:make-energy-function"),true,&chem::EnergyFunction_O::make,"(matter force_field &optional active-atoms)");
}
NOINLINE void expose_function_1075_helper() {
  expose_function(core::magic_name("chem__map_atoms"),true,&chem::chem__map_atoms,"");
}
NOINLINE void expose_function_1076_helper() {
  expose_function(core::magic_name("chem__map_residues"),true,&chem::chem__map_residues,"");
}
NOINLINE void expose_function_1077_helper() {
  expose_function(core::magic_name("chem__map_molecules"),true,&chem::chem__map_molecules,"");
}
NOINLINE void expose_function_1078_helper() {
  expose_function(core::magic_name("chem__readPdbMonomerConnectivityDatabase"),true,&chem::chem__readPdbMonomerConnectivityDatabase,"");
}
NOINLINE void expose_function_1079_helper() {
  expose_function(core::magic_name("CHEM:make-minimizer"),true,&chem::Minimizer_O::make,"(&key matter force_field energy_function)");
}
NOINLINE void expose_function_1080_helper() {
  expose_function(core::magic_name("CHEM:elementFromAtomNameString"),true,&chem::elementFromAtomNameString,"");
}
NOINLINE void expose_function_1081_helper() {
  expose_function(core::magic_name("CHEM:elementFromAtomNameStringCaseInsensitive"),true,&chem::elementFromAtomNameStringCaseInsensitive,"");
}
NOINLINE void expose_function_1082_helper() {
  expose_function(core::magic_name("CHEM:vdwRadiusForElement"),true,&chem::vdwRadiusForElement,"");
}
NOINLINE void expose_function_1083_helper() {
  expose_function(core::magic_name("CHEM:make-atom"),true,&chem::Atom_O::make,"");
}
NOINLINE void expose_function_1084_helper() {
  expose_function(core::magic_name("CHEM:make-residue"),true,&chem::Residue_O::make,"");
}
NOINLINE void expose_function_1085_helper() {
  expose_function(core::magic_name("CHEM:make-molecule"),true,&chem::Molecule_O::make,"(&optional (name nil))");
}
NOINLINE void expose_function_1086_helper() {
  expose_function(core::magic_name("CHEM:make-aggregate"),true,&chem::Aggregate_O::make,"(&optional (name nil))");
}
NOINLINE void expose_function_1087_helper() {
  expose_function(core::magic_name("chem__make_fftypes_db"),true,&chem::chem__make_fftypes_db,"");
}
#endif // EXPOSE_FUNCTION_BINDINGS_HELPERS
#ifdef EXPOSE_FUNCTION_BINDINGS
  expose_function_1_helper();
  expose_function_2_helper();
  expose_function_3_helper();
  expose_function_4_helper();
  expose_function_5_helper();
  expose_function_6_helper();
  expose_function_7_helper();
  expose_function_8_helper();
  expose_function_9_helper();
  expose_function_10_helper();
  expose_function_11_helper();
  expose_function_12_helper();
  expose_function_13_helper();
  expose_function_14_helper();
  expose_function_15_helper();
  expose_function_16_helper();
  expose_function_17_helper();
  expose_function_18_helper();
  expose_function_19_helper();
  expose_function_20_helper();
  expose_function_21_helper();
  expose_function_22_helper();
  expose_function_23_helper();
  expose_function_24_helper();
  expose_function_25_helper();
  expose_function_26_helper();
  expose_function_27_helper();
  expose_function_28_helper();
  expose_function_29_helper();
  expose_function_30_helper();
  expose_function_31_helper();
  expose_function_32_helper();
  expose_function_33_helper();
  expose_function_34_helper();
  expose_function_35_helper();
  expose_function_36_helper();
  expose_function_37_helper();
  expose_function_38_helper();
  expose_function_39_helper();
  expose_function_40_helper();
  expose_function_41_helper();
  expose_function_42_helper();
  expose_function_43_helper();
  expose_function_44_helper();
  expose_function_45_helper();
  expose_function_46_helper();
  expose_function_47_helper();
  expose_function_48_helper();
  expose_function_49_helper();
  expose_function_50_helper();
  expose_function_51_helper();
  expose_function_52_helper();
  expose_function_53_helper();
  expose_function_54_helper();
  expose_function_55_helper();
  expose_function_56_helper();
  expose_function_57_helper();
  expose_function_58_helper();
  expose_function_59_helper();
  expose_function_60_helper();
  expose_function_61_helper();
  expose_function_62_helper();
  expose_function_63_helper();
  expose_function_64_helper();
  expose_function_65_helper();
  expose_function_66_helper();
  expose_function_67_helper();
  expose_function_68_helper();
  expose_function_69_helper();
  expose_function_70_helper();
  expose_function_71_helper();
  expose_function_72_helper();
  expose_function_73_helper();
  expose_function_74_helper();
  expose_function_75_helper();
  expose_function_76_helper();
  expose_function_77_helper();
  expose_function_78_helper();
  expose_function_79_helper();
  expose_function_80_helper();
  expose_function_81_helper();
  expose_function_82_helper();
  expose_function_83_helper();
  expose_function_84_helper();
  expose_function_85_helper();
  expose_function_86_helper();
  expose_function_87_helper();
  expose_function_88_helper();
  expose_function_89_helper();
  expose_function_90_helper();
  expose_function_91_helper();
  expose_function_92_helper();
  expose_function_93_helper();
  expose_function_94_helper();
  expose_function_95_helper();
  expose_function_96_helper();
  expose_function_97_helper();
  expose_function_98_helper();
  expose_function_99_helper();
  expose_function_100_helper();
  expose_function_101_helper();
  expose_function_102_helper();
  expose_function_103_helper();
  expose_function_104_helper();
  expose_function_105_helper();
  expose_function_106_helper();
  expose_function_107_helper();
  expose_function_108_helper();
  expose_function_109_helper();
  expose_function_110_helper();
  expose_function_111_helper();
  expose_function_112_helper();
  expose_function_113_helper();
  expose_function_114_helper();
  expose_function_115_helper();
  expose_function_116_helper();
  expose_function_117_helper();
  expose_function_118_helper();
  expose_function_119_helper();
  expose_function_120_helper();
  expose_function_121_helper();
  expose_function_122_helper();
  expose_function_123_helper();
  expose_function_124_helper();
  expose_function_125_helper();
  expose_function_126_helper();
  expose_function_127_helper();
  expose_function_128_helper();
  expose_function_129_helper();
  expose_function_130_helper();
  expose_function_131_helper();
  expose_function_132_helper();
  expose_function_133_helper();
  expose_function_134_helper();
  expose_function_135_helper();
  expose_function_136_helper();
  expose_function_137_helper();
  expose_function_138_helper();
  expose_function_139_helper();
  expose_function_140_helper();
  expose_function_141_helper();
  expose_function_142_helper();
  expose_function_143_helper();
  expose_function_144_helper();
  expose_function_145_helper();
  expose_function_146_helper();
  expose_function_147_helper();
  expose_function_148_helper();
  expose_function_149_helper();
  expose_function_150_helper();
  expose_function_151_helper();
  expose_function_152_helper();
  expose_function_153_helper();
  expose_function_154_helper();
  expose_function_155_helper();
  expose_function_156_helper();
  expose_function_157_helper();
  expose_function_158_helper();
  expose_function_159_helper();
  expose_function_160_helper();
  expose_function_161_helper();
  expose_function_162_helper();
  expose_function_163_helper();
  expose_function_164_helper();
  expose_function_165_helper();
  expose_function_166_helper();
  expose_function_167_helper();
  expose_function_168_helper();
  expose_function_169_helper();
  expose_function_170_helper();
  expose_function_171_helper();
  expose_function_172_helper();
  expose_function_173_helper();
  expose_function_174_helper();
  expose_function_175_helper();
  expose_function_176_helper();
  expose_function_177_helper();
  expose_function_178_helper();
  expose_function_179_helper();
  expose_function_180_helper();
  expose_function_181_helper();
  expose_function_182_helper();
  expose_function_183_helper();
  expose_function_184_helper();
  expose_function_185_helper();
  expose_function_186_helper();
  expose_function_187_helper();
  expose_function_188_helper();
  expose_function_189_helper();
  expose_function_190_helper();
  expose_function_191_helper();
  expose_function_192_helper();
  expose_function_193_helper();
  expose_function_194_helper();
  expose_function_195_helper();
  expose_function_196_helper();
  expose_function_197_helper();
  expose_function_198_helper();
  expose_function_199_helper();
  expose_function_200_helper();
  expose_function_201_helper();
  expose_function_202_helper();
  expose_function_203_helper();
  expose_function_204_helper();
  expose_function_205_helper();
  expose_function_206_helper();
  expose_function_207_helper();
  expose_function_208_helper();
  expose_function_209_helper();
  expose_function_210_helper();
  expose_function_211_helper();
  expose_function_212_helper();
  expose_function_213_helper();
  expose_function_214_helper();
  expose_function_215_helper();
  expose_function_216_helper();
  expose_function_217_helper();
  expose_function_218_helper();
  expose_function_219_helper();
  expose_function_220_helper();
  expose_function_221_helper();
  expose_function_222_helper();
  expose_function_223_helper();
  expose_function_224_helper();
  expose_function_225_helper();
  expose_function_226_helper();
  expose_function_227_helper();
  expose_function_228_helper();
  expose_function_229_helper();
  expose_function_230_helper();
  expose_function_231_helper();
  expose_function_232_helper();
  expose_function_233_helper();
  expose_function_234_helper();
  expose_function_235_helper();
  expose_function_236_helper();
  expose_function_237_helper();
  expose_function_238_helper();
  expose_function_239_helper();
  expose_function_240_helper();
  expose_function_241_helper();
  expose_function_242_helper();
  expose_function_243_helper();
  expose_function_244_helper();
  expose_function_245_helper();
  expose_function_246_helper();
  expose_function_247_helper();
  expose_function_248_helper();
  expose_function_249_helper();
  expose_function_250_helper();
  expose_function_251_helper();
  expose_function_252_helper();
  expose_function_253_helper();
  expose_function_254_helper();
  expose_function_255_helper();
  expose_function_256_helper();
  expose_function_257_helper();
  expose_function_258_helper();
  expose_function_259_helper();
  expose_function_260_helper();
  expose_function_261_helper();
  expose_function_262_helper();
  expose_function_263_helper();
  expose_function_264_helper();
  expose_function_265_helper();
  expose_function_266_helper();
  expose_function_267_helper();
  expose_function_268_helper();
  expose_function_269_helper();
  expose_function_270_helper();
  expose_function_271_helper();
  expose_function_272_helper();
  expose_function_273_helper();
  expose_function_274_helper();
  expose_function_275_helper();
  expose_function_276_helper();
  expose_function_277_helper();
  expose_function_278_helper();
  expose_function_279_helper();
  expose_function_280_helper();
  expose_function_281_helper();
  expose_function_282_helper();
  expose_function_283_helper();
  expose_function_284_helper();
  expose_function_285_helper();
  expose_function_286_helper();
  expose_function_287_helper();
  expose_function_288_helper();
  expose_function_289_helper();
  expose_function_290_helper();
  expose_function_291_helper();
  expose_function_292_helper();
  expose_function_293_helper();
  expose_function_294_helper();
  expose_function_295_helper();
  expose_function_296_helper();
  expose_function_297_helper();
  expose_function_298_helper();
  expose_function_299_helper();
  expose_function_300_helper();
  expose_function_301_helper();
  expose_function_302_helper();
  expose_function_303_helper();
  expose_function_304_helper();
  expose_function_305_helper();
  expose_function_306_helper();
  expose_function_307_helper();
  expose_function_308_helper();
  expose_function_309_helper();
  expose_function_310_helper();
  expose_function_311_helper();
  expose_function_312_helper();
  expose_function_313_helper();
  expose_function_314_helper();
  expose_function_315_helper();
  expose_function_316_helper();
  expose_function_317_helper();
  expose_function_318_helper();
  expose_function_319_helper();
  expose_function_320_helper();
  expose_function_321_helper();
  expose_function_322_helper();
  expose_function_323_helper();
  expose_function_324_helper();
  expose_function_325_helper();
  expose_function_326_helper();
  expose_function_327_helper();
  expose_function_328_helper();
  expose_function_329_helper();
  expose_function_330_helper();
  expose_function_331_helper();
  expose_function_332_helper();
  expose_function_333_helper();
  expose_function_334_helper();
  expose_function_335_helper();
  expose_function_336_helper();
  expose_function_337_helper();
  expose_function_338_helper();
  expose_function_339_helper();
  expose_function_340_helper();
  expose_function_341_helper();
  expose_function_342_helper();
  expose_function_343_helper();
  expose_function_344_helper();
  expose_function_345_helper();
  expose_function_346_helper();
  expose_function_347_helper();
  expose_function_348_helper();
  expose_function_349_helper();
  expose_function_350_helper();
  expose_function_351_helper();
  expose_function_352_helper();
  expose_function_353_helper();
  expose_function_354_helper();
  expose_function_355_helper();
  expose_function_356_helper();
  expose_function_357_helper();
  expose_function_358_helper();
  expose_function_359_helper();
  expose_function_360_helper();
  expose_function_361_helper();
  expose_function_362_helper();
  expose_function_363_helper();
  expose_function_364_helper();
  expose_function_365_helper();
  expose_function_366_helper();
  expose_function_367_helper();
  expose_function_368_helper();
  expose_function_369_helper();
  expose_function_370_helper();
  expose_function_371_helper();
  expose_function_372_helper();
  expose_function_373_helper();
  expose_function_374_helper();
  expose_function_375_helper();
  expose_function_376_helper();
  expose_function_377_helper();
  expose_function_378_helper();
  expose_function_379_helper();
  expose_function_380_helper();
  expose_function_381_helper();
  expose_function_382_helper();
  expose_function_383_helper();
  expose_function_384_helper();
  expose_function_385_helper();
  expose_function_386_helper();
  expose_function_387_helper();
  expose_function_388_helper();
  expose_function_389_helper();
  expose_function_390_helper();
  expose_function_391_helper();
  expose_function_392_helper();
  expose_function_393_helper();
  expose_function_394_helper();
  expose_function_395_helper();
  expose_function_396_helper();
  expose_function_397_helper();
  expose_function_398_helper();
  expose_function_399_helper();
  expose_function_400_helper();
  expose_function_401_helper();
  expose_function_402_helper();
  expose_function_403_helper();
  expose_function_404_helper();
  expose_function_405_helper();
  expose_function_406_helper();
  expose_function_407_helper();
  expose_function_408_helper();
  expose_function_409_helper();
  expose_function_410_helper();
  expose_function_411_helper();
  expose_function_412_helper();
  expose_function_413_helper();
  expose_function_414_helper();
  expose_function_415_helper();
  expose_function_416_helper();
  expose_function_417_helper();
  expose_function_418_helper();
  expose_function_419_helper();
  expose_function_420_helper();
  expose_function_421_helper();
  expose_function_422_helper();
  expose_function_423_helper();
  expose_function_424_helper();
  expose_function_425_helper();
  expose_function_426_helper();
  expose_function_427_helper();
  expose_function_428_helper();
  expose_function_429_helper();
  expose_function_430_helper();
  expose_function_431_helper();
  expose_function_432_helper();
  expose_function_433_helper();
  expose_function_434_helper();
  expose_function_435_helper();
  expose_function_436_helper();
  expose_function_437_helper();
  expose_function_438_helper();
  expose_function_439_helper();
  expose_function_440_helper();
  expose_function_441_helper();
  expose_function_442_helper();
  expose_function_443_helper();
  expose_function_444_helper();
  expose_function_445_helper();
  expose_function_446_helper();
  expose_function_447_helper();
  expose_function_448_helper();
  expose_function_449_helper();
  expose_function_450_helper();
  expose_function_451_helper();
  expose_function_452_helper();
  expose_function_453_helper();
  expose_function_454_helper();
  expose_function_455_helper();
  expose_function_456_helper();
  expose_function_457_helper();
  expose_function_458_helper();
  expose_function_459_helper();
  expose_function_460_helper();
  expose_function_461_helper();
  expose_function_462_helper();
  expose_function_463_helper();
  expose_function_464_helper();
  expose_function_465_helper();
  expose_function_466_helper();
  expose_function_467_helper();
  expose_function_468_helper();
  expose_function_469_helper();
  expose_function_470_helper();
  expose_function_471_helper();
  expose_function_472_helper();
  expose_function_473_helper();
  expose_function_474_helper();
  expose_function_475_helper();
  expose_function_476_helper();
  expose_function_477_helper();
  expose_function_478_helper();
  expose_function_479_helper();
  expose_function_480_helper();
  expose_function_481_helper();
  expose_function_482_helper();
  expose_function_483_helper();
  expose_function_484_helper();
  expose_function_485_helper();
  expose_function_486_helper();
  expose_function_487_helper();
  expose_function_488_helper();
  expose_function_489_helper();
  expose_function_490_helper();
  expose_function_491_helper();
  expose_function_492_helper();
  expose_function_493_helper();
  expose_function_494_helper();
  expose_function_495_helper();
  expose_function_496_helper();
  expose_function_497_helper();
  expose_function_498_helper();
  expose_function_499_helper();
  expose_function_500_helper();
  expose_function_501_helper();
  expose_function_502_helper();
  expose_function_503_helper();
  expose_function_504_helper();
  expose_function_505_helper();
  expose_function_506_helper();
  expose_function_507_helper();
  expose_function_508_helper();
  expose_function_509_helper();
  expose_function_510_helper();
  expose_function_511_helper();
  expose_function_512_helper();
  expose_function_513_helper();
  expose_function_514_helper();
  expose_function_515_helper();
  expose_function_516_helper();
  expose_function_517_helper();
  expose_function_518_helper();
  expose_function_519_helper();
  expose_function_520_helper();
  expose_function_521_helper();
  expose_function_522_helper();
  expose_function_523_helper();
  expose_function_524_helper();
  expose_function_525_helper();
  expose_function_526_helper();
  expose_function_527_helper();
  expose_function_528_helper();
  expose_function_529_helper();
  expose_function_530_helper();
  expose_function_531_helper();
  expose_function_532_helper();
  expose_function_533_helper();
  expose_function_534_helper();
  expose_function_535_helper();
  expose_function_536_helper();
  expose_function_537_helper();
  expose_function_538_helper();
  expose_function_539_helper();
  expose_function_540_helper();
  expose_function_541_helper();
  expose_function_542_helper();
  expose_function_543_helper();
  expose_function_544_helper();
  expose_function_545_helper();
  expose_function_546_helper();
  expose_function_547_helper();
  expose_function_548_helper();
  expose_function_549_helper();
  expose_function_550_helper();
  expose_function_551_helper();
  expose_function_552_helper();
  expose_function_553_helper();
  expose_function_554_helper();
  expose_function_555_helper();
  expose_function_556_helper();
  expose_function_557_helper();
  expose_function_558_helper();
  expose_function_559_helper();
  expose_function_560_helper();
  expose_function_561_helper();
  expose_function_562_helper();
  expose_function_563_helper();
  expose_function_564_helper();
  expose_function_565_helper();
  expose_function_566_helper();
  expose_function_567_helper();
  expose_function_568_helper();
  expose_function_569_helper();
  expose_function_570_helper();
  expose_function_571_helper();
  expose_function_572_helper();
  expose_function_573_helper();
  expose_function_574_helper();
  expose_function_575_helper();
  expose_function_576_helper();
  expose_function_577_helper();
  expose_function_578_helper();
  expose_function_579_helper();
  expose_function_580_helper();
  expose_function_581_helper();
  expose_function_582_helper();
  expose_function_583_helper();
  expose_function_584_helper();
  expose_function_585_helper();
  expose_function_586_helper();
  expose_function_587_helper();
  expose_function_588_helper();
  expose_function_589_helper();
  expose_function_590_helper();
  expose_function_591_helper();
  expose_function_592_helper();
  expose_function_593_helper();
  expose_function_594_helper();
  expose_function_595_helper();
  expose_function_596_helper();
  expose_function_597_helper();
  expose_function_598_helper();
  expose_function_599_helper();
  expose_function_600_helper();
  expose_function_601_helper();
  expose_function_602_helper();
  expose_function_603_helper();
  expose_function_604_helper();
  expose_function_605_helper();
  expose_function_606_helper();
  expose_function_607_helper();
  expose_function_608_helper();
  expose_function_609_helper();
  expose_function_610_helper();
  expose_function_611_helper();
  expose_function_612_helper();
  expose_function_613_helper();
  expose_function_614_helper();
  expose_function_615_helper();
  expose_function_616_helper();
  expose_function_617_helper();
  expose_function_618_helper();
  expose_function_619_helper();
  expose_function_620_helper();
  expose_function_621_helper();
  expose_function_622_helper();
  expose_function_623_helper();
  expose_function_624_helper();
  expose_function_625_helper();
  expose_function_626_helper();
  expose_function_627_helper();
  expose_function_628_helper();
  expose_function_629_helper();
  expose_function_630_helper();
  expose_function_631_helper();
  expose_function_632_helper();
  expose_function_633_helper();
  expose_function_634_helper();
  expose_function_635_helper();
  expose_function_636_helper();
  expose_function_637_helper();
  expose_function_638_helper();
  expose_function_639_helper();
  expose_function_640_helper();
  expose_function_641_helper();
  expose_function_642_helper();
  expose_function_643_helper();
  expose_function_644_helper();
  expose_function_645_helper();
  expose_function_646_helper();
  expose_function_647_helper();
  expose_function_648_helper();
  expose_function_649_helper();
  expose_function_650_helper();
  expose_function_651_helper();
  expose_function_652_helper();
  expose_function_653_helper();
  expose_function_654_helper();
  expose_function_655_helper();
  expose_function_656_helper();
  expose_function_657_helper();
  expose_function_658_helper();
  expose_function_659_helper();
  expose_function_660_helper();
  expose_function_661_helper();
  expose_function_662_helper();
  expose_function_663_helper();
  expose_function_664_helper();
  expose_function_665_helper();
  expose_function_666_helper();
  expose_function_667_helper();
  expose_function_668_helper();
  expose_function_669_helper();
  expose_function_670_helper();
  expose_function_671_helper();
  expose_function_672_helper();
  expose_function_673_helper();
  expose_function_674_helper();
  expose_function_675_helper();
  expose_function_676_helper();
  expose_function_677_helper();
  expose_function_678_helper();
  expose_function_679_helper();
  expose_function_680_helper();
  expose_function_681_helper();
  expose_function_682_helper();
  expose_function_683_helper();
  expose_function_684_helper();
  expose_function_685_helper();
  expose_function_686_helper();
  expose_function_687_helper();
  expose_function_688_helper();
  expose_function_689_helper();
  expose_function_690_helper();
  expose_function_691_helper();
  expose_function_692_helper();
  expose_function_693_helper();
  expose_function_694_helper();
  expose_function_695_helper();
  expose_function_696_helper();
  expose_function_697_helper();
  expose_function_698_helper();
  expose_function_699_helper();
  expose_function_700_helper();
  expose_function_701_helper();
  expose_function_702_helper();
  expose_function_703_helper();
  expose_function_704_helper();
  expose_function_705_helper();
  expose_function_706_helper();
  expose_function_707_helper();
  expose_function_708_helper();
  expose_function_709_helper();
  expose_function_710_helper();
  expose_function_711_helper();
  expose_function_712_helper();
  expose_function_713_helper();
  expose_function_714_helper();
  expose_function_715_helper();
  expose_function_716_helper();
  expose_function_717_helper();
  expose_function_718_helper();
  expose_function_719_helper();
  expose_function_720_helper();
  expose_function_721_helper();
  expose_function_722_helper();
  expose_function_723_helper();
  expose_function_724_helper();
  expose_function_725_helper();
  expose_function_726_helper();
  expose_function_727_helper();
  expose_function_728_helper();
  expose_function_729_helper();
  expose_function_730_helper();
  expose_function_731_helper();
  expose_function_732_helper();
  expose_function_733_helper();
  expose_function_734_helper();
  expose_function_735_helper();
  expose_function_736_helper();
  expose_function_737_helper();
  expose_function_738_helper();
  expose_function_739_helper();
  expose_function_740_helper();
  expose_function_741_helper();
  expose_function_742_helper();
  expose_function_743_helper();
  expose_function_744_helper();
  expose_function_745_helper();
  expose_function_746_helper();
  expose_function_747_helper();
  expose_function_748_helper();
  expose_function_749_helper();
  expose_function_750_helper();
  expose_function_751_helper();
  expose_function_752_helper();
  expose_function_753_helper();
  expose_function_754_helper();
  expose_function_755_helper();
  expose_function_756_helper();
  expose_function_757_helper();
  expose_function_758_helper();
  expose_function_759_helper();
  expose_function_760_helper();
  expose_function_761_helper();
  expose_function_762_helper();
  expose_function_763_helper();
  expose_function_764_helper();
  expose_function_765_helper();
  expose_function_766_helper();
  expose_function_767_helper();
  expose_function_768_helper();
  expose_function_769_helper();
  expose_function_770_helper();
  expose_function_771_helper();
  expose_function_772_helper();
  expose_function_773_helper();
  expose_function_774_helper();
  expose_function_775_helper();
  expose_function_776_helper();
  expose_function_777_helper();
  expose_function_778_helper();
  expose_function_779_helper();
  expose_function_780_helper();
  expose_function_781_helper();
  expose_function_782_helper();
  expose_function_783_helper();
  expose_function_784_helper();
  expose_function_785_helper();
  expose_function_786_helper();
  expose_function_787_helper();
  expose_function_788_helper();
  expose_function_789_helper();
  expose_function_790_helper();
  expose_function_791_helper();
  expose_function_792_helper();
  expose_function_793_helper();
  expose_function_794_helper();
  expose_function_795_helper();
  expose_function_796_helper();
  expose_function_797_helper();
  expose_function_798_helper();
  expose_function_799_helper();
  expose_function_800_helper();
  expose_function_801_helper();
  expose_function_802_helper();
  expose_function_803_helper();
  expose_function_804_helper();
  expose_function_805_helper();
  expose_function_806_helper();
  expose_function_807_helper();
  expose_function_808_helper();
  expose_function_809_helper();
  expose_function_810_helper();
  expose_function_811_helper();
  expose_function_812_helper();
  expose_function_813_helper();
  expose_function_814_helper();
  expose_function_815_helper();
  expose_function_816_helper();
  expose_function_817_helper();
  expose_function_818_helper();
  expose_function_819_helper();
  expose_function_820_helper();
  expose_function_821_helper();
  expose_function_822_helper();
  expose_function_823_helper();
  expose_function_824_helper();
  expose_function_825_helper();
  expose_function_826_helper();
  expose_function_827_helper();
  expose_function_828_helper();
  expose_function_829_helper();
  expose_function_830_helper();
  expose_function_831_helper();
  expose_function_832_helper();
  expose_function_833_helper();
  expose_function_834_helper();
  expose_function_835_helper();
  expose_function_836_helper();
  expose_function_837_helper();
  expose_function_838_helper();
  expose_function_839_helper();
  expose_function_840_helper();
  expose_function_841_helper();
  expose_function_842_helper();
  expose_function_843_helper();
  expose_function_844_helper();
  expose_function_845_helper();
  expose_function_846_helper();
  expose_function_847_helper();
  expose_function_848_helper();
  expose_function_849_helper();
  expose_function_850_helper();
  expose_function_851_helper();
  expose_function_852_helper();
  expose_function_853_helper();
  expose_function_854_helper();
  expose_function_855_helper();
  expose_function_856_helper();
  expose_function_857_helper();
  expose_function_858_helper();
  expose_function_859_helper();
  expose_function_860_helper();
  expose_function_861_helper();
  expose_function_862_helper();
  expose_function_863_helper();
  expose_function_864_helper();
  expose_function_865_helper();
  expose_function_866_helper();
  expose_function_867_helper();
  expose_function_868_helper();
  expose_function_869_helper();
  expose_function_870_helper();
  expose_function_871_helper();
  expose_function_872_helper();
  expose_function_873_helper();
  expose_function_874_helper();
  expose_function_875_helper();
  expose_function_876_helper();
  expose_function_877_helper();
  expose_function_878_helper();
  expose_function_879_helper();
  expose_function_880_helper();
  expose_function_881_helper();
  expose_function_882_helper();
  expose_function_883_helper();
  expose_function_884_helper();
  expose_function_885_helper();
  expose_function_886_helper();
  expose_function_887_helper();
  expose_function_888_helper();
  expose_function_889_helper();
  expose_function_890_helper();
  expose_function_891_helper();
  expose_function_892_helper();
  expose_function_893_helper();
  expose_function_894_helper();
  expose_function_895_helper();
  expose_function_896_helper();
  expose_function_897_helper();
  expose_function_898_helper();
  expose_function_899_helper();
  expose_function_900_helper();
  expose_function_901_helper();
  expose_function_902_helper();
  expose_function_903_helper();
  expose_function_904_helper();
  expose_function_905_helper();
  expose_function_906_helper();
  expose_function_907_helper();
  expose_function_908_helper();
  expose_function_909_helper();
  expose_function_910_helper();
  expose_function_911_helper();
  expose_function_912_helper();
  expose_function_913_helper();
  expose_function_914_helper();
  expose_function_915_helper();
  expose_function_916_helper();
  expose_function_917_helper();
  expose_function_918_helper();
  expose_function_919_helper();
  expose_function_920_helper();
  expose_function_921_helper();
  expose_function_922_helper();
  expose_function_923_helper();
  expose_function_924_helper();
  expose_function_925_helper();
  expose_function_926_helper();
  expose_function_927_helper();
  expose_function_928_helper();
  expose_function_929_helper();
  expose_function_930_helper();
  expose_function_931_helper();
  expose_function_932_helper();
  expose_function_933_helper();
  expose_function_934_helper();
  expose_function_935_helper();
  expose_function_936_helper();
  expose_function_937_helper();
  expose_function_938_helper();
  expose_function_939_helper();
  expose_function_940_helper();
  expose_function_941_helper();
  expose_function_942_helper();
  expose_function_943_helper();
  expose_function_944_helper();
  expose_function_945_helper();
  expose_function_946_helper();
  expose_function_947_helper();
  expose_function_948_helper();
  expose_function_949_helper();
  expose_function_950_helper();
  expose_function_951_helper();
  expose_function_952_helper();
  expose_function_953_helper();
  expose_function_954_helper();
  expose_function_955_helper();
  expose_function_956_helper();
  expose_function_957_helper();
  expose_function_958_helper();
  expose_function_959_helper();
  expose_function_960_helper();
  expose_function_961_helper();
  expose_function_962_helper();
  expose_function_963_helper();
  expose_function_964_helper();
  expose_function_965_helper();
  expose_function_966_helper();
  expose_function_967_helper();
  expose_function_968_helper();
  expose_function_969_helper();
  expose_function_970_helper();
  expose_function_971_helper();
  expose_function_972_helper();
  expose_function_973_helper();
  expose_function_974_helper();
  expose_function_975_helper();
  expose_function_976_helper();
  expose_function_977_helper();
  expose_function_978_helper();
  expose_function_979_helper();
  expose_function_980_helper();
  expose_function_981_helper();
  expose_function_982_helper();
  expose_function_983_helper();
  expose_function_984_helper();
  expose_function_985_helper();
  expose_function_986_helper();
  expose_function_987_helper();
  expose_function_988_helper();
  expose_function_989_helper();
  expose_function_990_helper();
  expose_function_991_helper();
  expose_function_992_helper();
  expose_function_993_helper();
  expose_function_994_helper();
  expose_function_995_helper();
  expose_function_996_helper();
  expose_function_997_helper();
  expose_function_998_helper();
  expose_function_999_helper();
  expose_function_1000_helper();
  expose_function_1001_helper();
  expose_function_1002_helper();
  expose_function_1003_helper();
  expose_function_1004_helper();
  expose_function_1005_helper();
  expose_function_1006_helper();
  expose_function_1007_helper();
  expose_function_1008_helper();
  expose_function_1009_helper();
  expose_function_1010_helper();
  expose_function_1011_helper();
  expose_function_1012_helper();
  expose_function_1013_helper();
  expose_function_1014_helper();
  expose_function_1015_helper();
  expose_function_1016_helper();
  expose_function_1017_helper();
  expose_function_1018_helper();
  expose_function_1019_helper();
  expose_function_1020_helper();
  expose_function_1021_helper();
  expose_function_1022_helper();
  expose_function_1023_helper();
  expose_function_1024_helper();
  expose_function_1025_helper();
  expose_function_1026_helper();
  expose_function_1027_helper();
  expose_function_1028_helper();
  expose_function_1029_helper();
  expose_function_1030_helper();
  expose_function_1031_helper();
  expose_function_1032_helper();
  expose_function_1033_helper();
  expose_function_1034_helper();
  expose_function_1035_helper();
  expose_function_1036_helper();
  expose_function_1037_helper();
  expose_function_1038_helper();
  expose_function_1039_helper();
  expose_function_1040_helper();
  expose_function_1041_helper();
  expose_function_1042_helper();
  expose_function_1043_helper();
  expose_function_1044_helper();
  expose_function_1045_helper();
  expose_function_1046_helper();
  expose_function_1047_helper();
  expose_function_1048_helper();
  expose_function_1049_helper();
  expose_function_1050_helper();
  expose_function_1051_helper();
  expose_function_1052_helper();
  expose_function_1053_helper();
  expose_function_1054_helper();
  expose_function_1055_helper();
  expose_function_1056_helper();
  expose_function_1057_helper();
  expose_function_1058_helper();
  expose_function_1059_helper();
  expose_function_1060_helper();
  expose_function_1061_helper();
  expose_function_1062_helper();
  expose_function_1063_helper();
  expose_function_1064_helper();
  expose_function_1065_helper();
  expose_function_1066_helper();
  expose_function_1067_helper();
  expose_function_1068_helper();
  expose_function_1069_helper();
  expose_function_1070_helper();
  expose_function_1071_helper();
  expose_function_1072_helper();
  expose_function_1073_helper();
  expose_function_1074_helper();
  expose_function_1075_helper();
  expose_function_1076_helper();
  expose_function_1077_helper();
  expose_function_1078_helper();
  expose_function_1079_helper();
  expose_function_1080_helper();
  expose_function_1081_helper();
  expose_function_1082_helper();
  expose_function_1083_helper();
  expose_function_1084_helper();
  expose_function_1085_helper();
  expose_function_1086_helper();
  expose_function_1087_helper();
#endif // EXPOSE_FUNCTION_BINDINGS
