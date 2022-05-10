
if(NOT "D:/cpp/Sci/Flow/build/_deps/fmtlib-subbuild/fmtlib-populate-prefix/src/fmtlib-populate-stamp/fmtlib-populate-gitinfo.txt" IS_NEWER_THAN "D:/cpp/Sci/Flow/build/_deps/fmtlib-subbuild/fmtlib-populate-prefix/src/fmtlib-populate-stamp/fmtlib-populate-gitclone-lastrun.txt")
  message(STATUS "Avoiding repeated git clone, stamp file is up to date: 'D:/cpp/Sci/Flow/build/_deps/fmtlib-subbuild/fmtlib-populate-prefix/src/fmtlib-populate-stamp/fmtlib-populate-gitclone-lastrun.txt'")
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E rm -rf "D:/cpp/Sci/Flow/build/_deps/fmtlib-src"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: 'D:/cpp/Sci/Flow/build/_deps/fmtlib-src'")
endif()

# try the clone 3 times in case there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "D:/Program Files/Git/cmd/git.exe"  clone --no-checkout --config "advice.detachedHead=false" "https://github.com/fmtlib/fmt.git" "fmtlib-src"
    WORKING_DIRECTORY "D:/cpp/Sci/Flow/build/_deps"
    RESULT_VARIABLE error_code
    )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once:
          ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: 'https://github.com/fmtlib/fmt.git'")
endif()

execute_process(
  COMMAND "D:/Program Files/Git/cmd/git.exe"  checkout 5.3.0 --
  WORKING_DIRECTORY "D:/cpp/Sci/Flow/build/_deps/fmtlib-src"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: '5.3.0'")
endif()

set(init_submodules TRUE)
if(init_submodules)
  execute_process(
    COMMAND "D:/Program Files/Git/cmd/git.exe"  submodule update --recursive --init 
    WORKING_DIRECTORY "D:/cpp/Sci/Flow/build/_deps/fmtlib-src"
    RESULT_VARIABLE error_code
    )
endif()
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: 'D:/cpp/Sci/Flow/build/_deps/fmtlib-src'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy
    "D:/cpp/Sci/Flow/build/_deps/fmtlib-subbuild/fmtlib-populate-prefix/src/fmtlib-populate-stamp/fmtlib-populate-gitinfo.txt"
    "D:/cpp/Sci/Flow/build/_deps/fmtlib-subbuild/fmtlib-populate-prefix/src/fmtlib-populate-stamp/fmtlib-populate-gitclone-lastrun.txt"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: 'D:/cpp/Sci/Flow/build/_deps/fmtlib-subbuild/fmtlib-populate-prefix/src/fmtlib-populate-stamp/fmtlib-populate-gitclone-lastrun.txt'")
endif()

