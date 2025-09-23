<template>
  <span
    class="badge"
    :style="{ backgroundColor: color }"
  >
    <a :href="`/example_overview#${text}`">{{ text }}</a> 
  </span>
  &nbsp;
</template>

<script setup>
defineProps({
  color: {
    type: String,
    required: true
  },
  text: {
    type: String,
    required: true
  }
});
</script>

<style scoped>
.badge {
  display: inline-block;
  padding: 0.1em 0.25em;
  border-radius: 10%;
  font-size: 1em;
  user-select: none;
}
</style>

<!-- Example usage -->
<!--
<ExampleTag color="#3498db" text="Info" />
<ExampleTag color="green" text="Success" />
-->
