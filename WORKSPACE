
bind(
  name = "google_guava",
  actual = "@google_guava_19//jar",
)

maven_jar(
  name = "google_guava_19",
  artifact = "com.google.guava:guava:19.0",
)

bind(
  name = "google_java_common_geometry",
  actual = "@google_s2_geometry_java//:lib",
)

new_git_repository(
  name = "google_s2_geometry_java",
  remote = "https://github.com/google/s2-geometry-library-java.git",
  build_file = "BUILD.google_s2_geometry_java",
  commit = "c28f287b996c0cedc5516a0426fbd49f6c9611ec",
)

maven_jar(
  name = "junit4",
  artifact = "junit:junit:4.12",
)

maven_jar(
  name = "junit3",
  artifact = "junit:junit:3.8.1",
)

bind(
  name = "javax_nullable",
  actual = "@findbugs_jsr305//jar",
)

maven_jar(
  name = "findbugs_jsr305",
  artifact = "com.google.code.findbugs:jsr305:3.0.1",
)

