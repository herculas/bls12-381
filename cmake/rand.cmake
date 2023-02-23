INCLUDE(FetchContent)

FETCHCONTENT_DECLARE(
        Rand
        GIT_REPOSITORY https://github.com/herculas/Rand
        GIT_TAG main
)

FETCHCONTENT_MAKEAVAILABLE(Rand)